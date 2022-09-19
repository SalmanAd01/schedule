#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include "scheduler.h"

#define N_CORES     4
#define START_TEMP  2000.0
#define T_AMBIENT   0.00001
#define K           0.054

int hyperperiod;
int cschedlen;
static task_info_t * task_info;
static int n_tasks;
static int n_procs;

static schedule_list_element_t ** schedule;
static schedule_list_element_t ** neighbour;

static double sim_time = 0;
static long int sim_iterations = 0;

/* Helper functions:
 * double random_number(void);  generates random number in range [0, INT_MAX)
 * int gcd(int a, int b);       returns GCD of integers a and b
 * int lcm(int a, int b);       returns LCM of integers a and b
 * int LCM(int xs[], int l);    returns LCM of l integers xs
 */

int gcd(int a, int b)
{
    int t;
    while (b != 0)
    {
        a %= b;
        t = a;
        a = b;
        b = t;
    }
    return a;
}

int lcm(int a, int b)
{
    return ((a * b) / gcd(a, b));
}

int LCM(int xs[], int len)
{
    int i, l = 1;
    for (i = 0; i < len; i++)
    {
        l = lcm(l, xs[i]);
    }
    return l;
}


/* Task information related functions:
 * void init_task_info(void);   initializes task_info list
 */

void init_task_info()
{
    task_info = (task_info_t *)malloc(sizeof(task_info_t) * n_tasks);
    int i;
    for (i = 0; i < n_tasks; i++)
    {
        task_info[i].pid = i;
        scanf("%d", &task_info[i].period);
        scanf("%d", &task_info[i].deadline_rel);
        scanf("%d", &task_info[i].wcet);
        scanf("%d", &task_info[i].activation_time);
        scanf("%d", &task_info[i].energy_reqmt);
    }
}


/* Schedule and SA related functions:
 * void init_schedule(void);                            initializes schedule and neighbour data structures
 * double reduce_temp(double);                          reduces SA temperature
 * double cost_deadlines(schedule_list_element_t ** s); deadline heuristic
 * double cost_validity(schedule_list_element_t ** s);  validity heuristic
 * double cost_energy(schedule_list_element_t ** s);    energy consumption heuristic
 * double sigmoid(x);                                   returns sigmoid(x)
 * double util(s, n);                                   utility function to compare two schedules
 * void neighbour_gen(void);                            generates neighbour of current schedule
 * void display_schedule(void);                         displays schedule data structure
 * void display_neighbour(void);                        displays neighbour data structure (for debugging purposes only)
 */

void init_schedule(void)
{
    int i, j;
    cschedlen = 0;
    for (i = 0; i < n_tasks; i++)
    {
        cschedlen += hyperperiod / task_info[i].period;
    }
    schedule = (schedule_list_element_t **)malloc(sizeof(schedule_list_element_t *) * n_procs);
    neighbour = (schedule_list_element_t **)malloc(sizeof(schedule_list_element_t *) * n_procs);
    for (i = 0; i < n_procs; i++)
    {
        schedule[i] = (schedule_list_element_t *)malloc(sizeof(schedule_list_element_t) * n_procs * cschedlen);
        neighbour[i] = (schedule_list_element_t *)malloc(sizeof(schedule_list_element_t) * n_procs * cschedlen);
        for (j = 0; j < n_tasks; j++)
        {
            schedule[i][j].task = NULL;
            neighbour[i][j].task = NULL;
        }
    }
    int x = -1, y = -1;
    for (i = 0; i < n_tasks; i++)
    {
        for (j = 0; j < hyperperiod/task_info[i].period; j++)
        {
            x = (x + 1) % n_procs;
            if (x == 0) y++;
            schedule[x][y].task = &task_info[i];
        }
    }
}

double reduce_temp(double T)
{
    // T(t) = Ta + (T(0) - Ta) * e^(-kt)
    // Ta = 0; the surrounding temperature
    T = T_AMBIENT + (START_TEMP - T_AMBIENT) * exp(-K * sim_time);
    sim_time += 0.001;
    sim_iterations += 1;
    return T;
}

double cost_deadlines(schedule_list_element_t ** s)
{
    int i, j;
    int core_time;
    double deviation[n_tasks], dev_schedule = 0;
    for (i = 0; i < n_tasks; i++)
    {
        deviation[i] = 0;
    }
    for (i = 0; i < n_procs; i++)
    {
        j = 0;
        core_time = 0;
        while (s[i][j].task != NULL)
        {
            deviation[s[i][j].task->pid] += \
            (core_time + s[i][j].task->wcet) \
            - (core_time - (core_time % s[i][j].task->period) \
            + s[i][j].task->deadline_rel);
            core_time += s[i][j].task->wcet;
            j++;
        }
    }
    for (i = 0; i < n_tasks; i++)
    {
        deviation[i] /= (double)(hyperperiod / task_info[i].period);
        dev_schedule += deviation[i];
    }
    dev_schedule /= (double)n_tasks;
}

double cost_validity(schedule_list_element_t ** s)
{
    double validity = 1;
    int i, j, flag = 0, h;
    int core_time[n_procs], task_inst_counter[n_tasks];
    for (i = 0; i < n_procs; i++)
    {
        core_time[i] = 0;       
    }
    for (i = 0; i < n_tasks; i++)
    {
        task_inst_counter[i] = 0;
    }
    i = 0;
    j = 0;
    while (j < n_procs * cschedlen)
    {
        if (s[i][j].task == NULL) flag++;
        else
        {
            h = (core_time[i] / s[i][j].task->period) - task_inst_counter[s[i][j].task->pid];
            if (h != 0)
            {
                validity -= 1 / (double)cschedlen;
            }
            task_inst_counter[s[i][j].task->pid] += 1;
            core_time[i] += s[i][j].task->wcet;
        }
        i = (i + 1) % n_procs;
        if (i == 0) 
        {
            if (flag == n_procs) break;
            j++;
            flag = 0;
        }
    }
    return validity;
}

double cost_energy(schedule_list_element_t ** s)
{
    double energy_reqmt_mean = 0, dev = 0;
    int i, j, core_energy_reqmt[n_procs];
    for (i = 0; i < n_procs; i++)
    {
        j = 0;
        core_energy_reqmt[i] = 0;
        while (s[i][j].task != NULL)
        {
            core_energy_reqmt[i] += s[i][j].task->energy_reqmt;
            j++;
        }
        energy_reqmt_mean += core_energy_reqmt[i];
    }
    energy_reqmt_mean /= (float)n_procs;
    for (i = 0; i < n_procs; i++)
    {
    	dev += pow(core_energy_reqmt[i] - energy_reqmt_mean, 2);
    }
    dev = sqrt(dev / (float)n_procs);
}

double sigmoid(double x, double T)
{
    return (1 / (1 + exp(-x / T)));
}

double util(schedule_list_element_t ** s, schedule_list_element_t ** n, double T)
{
    return (sigmoid(cost_deadlines(n) - cost_deadlines(s), T) \
    + sigmoid(cost_validity(s) - cost_validity(n), T) \
    + sigmoid(cost_energy(n) - cost_energy(s), T));
}

void neighbour_gen()
{
    int i, j, k;
    double x;
    int r1, r2, r3, r4, r5, r6, r_temp;
    int ts[n_procs];
    schedule_list_element_t sel_temp;
    for (i = 0; i < n_procs; i++)
    {
        k = 0;
        for (j = 0; j < n_procs * cschedlen; j++)
        {
            if (schedule[i][j].task != NULL) k++;
            neighbour[i][j].task = schedule[i][j].task;
        }
        ts[i] = k;
    }
    for (i = 0; i < n_procs; i++)
    {
        r1 = rand() / (double) INT_MAX * ts[i];
        r2 = rand() / (double) INT_MAX * ts[i];
        if (r1 == r2) continue;
        else if (r1 > r2)
        {
            r_temp = r1;
            r1 = r2;
            r2 = r_temp;
        }
        while (r1 < r2)
        {
            sel_temp = neighbour[i][r1];
            neighbour[i][r1] = neighbour[i][r2];
            neighbour[i][r2] = sel_temp;
            r1++;
            r2--;
        }
    }
    r1 = rand() / (double) INT_MAX * n_procs;
    r2 = rand() / (double) INT_MAX * n_procs;
    if (r1 != r2)
    {
        r3 = rand() / (double) INT_MAX * ts[r1];
        r5 = r3;
        r4 = rand() / (double) INT_MAX * ts[r2];
        r6 = r4;
        while (r3 < ts[r1] || r4 < ts[r2])
        {
            sel_temp = neighbour[r1][r3];
            neighbour[r1][r3] = neighbour[r2][r4];
            neighbour[r2][r4] = sel_temp;
            r3++;
            r4++;
        }
        r3 = ts[r1];
        r4 = ts[r2];
        ts[r1] = r5 + (r4 - r6);
        ts[r2] = r6 + (r3 - r5);
    }
}

void display_schedule()
{
    int i, j;
    printf("\n\n");
    for (i = 0; i < n_procs; i++)
    {
        j = 0;
        printf("Core %d\n", i);
        while (schedule[i][j].task != NULL)
        {
            printf("    Task: ");
            printf("%d\n", schedule[i][j].task->pid);
            //printf("      Period: ");
            //printf("%d\n", schedule[i][j].task->period);
            j++;
        }
        printf("\n");
    }
}

void display_neighbour()
{
    int i, j;
    printf("\n\n");
    for (i = 0; i < n_procs; i++)
    {
        j = 0;
        printf("Core %d\n", i);
        while (neighbour[i][j].task != NULL)
        {
            printf("    Task: ");
            printf("%d\n", neighbour[i][j].task->pid);
            //printf("      Period: ");
            //printf("%d\n", neighbour[i][j].task->period);
            j++;
        }
        printf("\n");
    }
}

int main(void)
{
    int i, j, t_start = time(NULL);
    srand(time(NULL));
    double u, r, T = START_TEMP, xval;
    n_tasks = 3;
    n_procs = N_CORES;
    init_task_info();
    int task_periods[n_tasks];
    for (i = 0; i < n_tasks; i++)
    {
        task_periods[i] = task_info[i].period;
    }
    printf("\nTasks allocated\n");
    hyperperiod = LCM(task_periods, n_tasks);
    init_schedule();
    printf("Schedule initialized\n");
    schedule_list_element_t ** temp;
    display_schedule();
    printf("Stats:\n");
    while(sim_iterations < 125000)
    {
        neighbour_gen();
        r = rand() / (double) INT_MAX;
        u = util(schedule, neighbour, T);
        xval = u;
        /*if (sim_iterations == 75000 || sim_iterations == 75001)
        {
            display_schedule();
            display_neighbour();
        }*/
        printf("%F %F %F %F %F %F %F %F ", T, r, cost_deadlines(schedule), cost_deadlines(neighbour), cost_validity(schedule), cost_validity(neighbour), u, xval);
        if (cost_deadlines(neighbour) < cost_deadlines(schedule) \
        && cost_validity(neighbour) > cost_validity(schedule) \
        && cost_energy(neighbour) < cost_energy(schedule))
        {
            printf("Deterministic move\n");
            temp = schedule;
            schedule = neighbour;
            neighbour = temp;
            xval = util(schedule, neighbour, T);
        }
        else if (u >= r)
        {
            printf("Probabilistic move\n");
            temp = schedule;
            schedule = neighbour;
            neighbour = temp;
        }
        else
        {
            printf("No move\n");
        }
        T = reduce_temp(T);
    }
    printf("\nSchedule generated\n");
    display_schedule();
    int t_end = time(NULL);
    printf("\nSchedule cost: %F\nNeighbour cost: %F\nSimulated time: %F(min)\nIterations: %li\nActual simulation time: %u\n", cost_deadlines(schedule), cost_deadlines(neighbour), sim_time, sim_iterations, t_end - t_start);
    return 0;
}
