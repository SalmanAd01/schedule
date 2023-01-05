#ifndef SCHEDULER_H
#define SCHEDULER_H

typedef struct {
	int pid;
	int period;
	int deadline_rel;
	int activation_time;
	int wcet;
	int energy_reqmt;
} task_info_t;

typedef struct {
	task_info_t * task;
} schedule_list_element_t;

#endif
