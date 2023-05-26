#include <stdio.h>
#include <mpi.h>
#include <pthread.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define TASK_LIST_SIZE 1000
#define NON_SHARING_TASKS 50
#define TASK_LISTS_NUMBER 5
#define TASK_DIFFICULTY_COEFFICIENT 1000000

//tags
#define ASKING_TAG 123
#define TASK_NUMBER_REPLY_TAG 1111
#define TASKS_SENDING_TAG 123456

int rankOfProc, numOfProcs;
pthread_t threads[2];
pthread_mutex_t mutex;
int* taskWeights;
int tasksDone;
int tasksLeft;
double globalResult = 0;
double result = 0;

double totalDisbalanceRatio = 0;

void executeTaskList(int *tasks){
    int i = 0;
    while(tasksLeft > 0){
        //fprintf(stderr, "worker works: %d, tasks left: %d\n", rankOfProc, tasksLeft);
        pthread_mutex_lock(&mutex);
        int weightOfTask = tasks[i];
        pthread_mutex_unlock(&mutex);

        for(int j = 0; j < weightOfTask; j++){
            globalResult += sin(i) / 2;
        }
        pthread_mutex_lock(&mutex);
        tasksDone++;
        tasksLeft--;
        pthread_mutex_unlock(&mutex);
        i++;
    }
}

void generateTaskWeights(int* taskWeights, int currentIteration){
    for(int i = 0; i < TASK_LIST_SIZE; ++i)
        taskWeights[i] = abs(50 - i % 100) * abs(rankOfProc - (currentIteration % numOfProcs)) * TASK_DIFFICULTY_COEFFICIENT;
        //taskWeights[i] = (rankOfProc + 1) * 10000000;
}

void askForTasks(int* tasks){
    for(int j = 0; j < numOfProcs; j++){
        int additionalTasks;
        if(j != rankOfProc && j < numOfProcs){
            //fprintf(stderr, "asking: %d, numOfProcs: %d\n", j, numOfProcs);
            MPI_Send(&rankOfProc, 1, MPI_INT, j, ASKING_TAG, MPI_COMM_WORLD);
            MPI_Recv(&additionalTasks, 1, MPI_INT, j, TASK_NUMBER_REPLY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //fprintf(stderr, "additionalTasks got in %d: %d\n", rankOfProc, additionalTasks);
            if(additionalTasks != 0){
                pthread_mutex_lock(&mutex);
                MPI_Recv(&(taskWeights[tasksLeft]), additionalTasks, MPI_INT, j, TASKS_SENDING_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                tasksLeft += additionalTasks;
                pthread_mutex_unlock(&mutex);
            }
        }
    }    
}

void* executorsWork(void*){
    //fprintf(stderr, "Executor started\n");
    taskWeights = (int*)malloc(TASK_LIST_SIZE * numOfProcs * sizeof(int));
    double workTimeStart, workTimeEnd, timeForWork, minTime, maxTime;
   
    for(int currentTaskList = 0; currentTaskList < TASK_LISTS_NUMBER; ++currentTaskList){
        generateTaskWeights(taskWeights, currentTaskList);
        tasksDone = 0;
        tasksLeft = TASK_LIST_SIZE;
        workTimeStart = MPI_Wtime();
        executeTaskList(taskWeights);
        for(int i = 0; i < 2; ++i){
            askForTasks(taskWeights);
            executeTaskList(taskWeights);
        }
        workTimeEnd = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        timeForWork = workTimeEnd - workTimeStart;
        MPI_Reduce(&timeForWork, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&timeForWork, &minTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&globalResult, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        //printf("Rank: %d, Tasks done: %d, time: %lf\n", rankOfProc, tasksDone, timeForWork);
        if(rankOfProc == 0){
            double disbalance = maxTime - minTime;
            double disbalanceRatio = (disbalance / maxTime) * 100;
            totalDisbalanceRatio += disbalanceRatio;
            fprintf(stderr, "Max time: %lf, min time: %lf, result: %lf Disbalance: %lf Disbalance ratio: %lf\n", maxTime, minTime, result, disbalance, disbalanceRatio);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rankOfProc == 0){
        fprintf(stderr, "Average disbalance ratio: %lf\n", totalDisbalanceRatio / TASK_LISTS_NUMBER);
    }
    //fprintf(stderr, "worker ended\n");
    int exitFlag = -1;
    MPI_Send(&exitFlag, 1, MPI_INT, rankOfProc, ASKING_TAG, MPI_COMM_WORLD);
    free(taskWeights);
    return NULL;
}

void* managersWork(void*){
    //fprintf(stderr, "Manager started\n");
    while(1){
        int rankOfClient;
        MPI_Recv(&rankOfClient, 1, MPI_INT, MPI_ANY_SOURCE, ASKING_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //fprintf(stderr, "received from: %d in %d\n", rankOfClient, rankOfProc);
        if(rankOfClient == -1){
            //fprintf(stderr, "manager ended\n");
            pthread_exit(EXIT_SUCCESS);
        }
        pthread_mutex_lock(&mutex);
        int tasksToShare;
        if(tasksLeft >= NON_SHARING_TASKS){
            tasksToShare = tasksLeft / 4;
            //fprintf(stderr, "Sending %d tasks from %d to %d\n", tasksToShare, rankOfProc, rankOfClient);
            MPI_Send(&tasksToShare, 1, MPI_INT, rankOfClient, TASK_NUMBER_REPLY_TAG, MPI_COMM_WORLD);
            MPI_Send(&(taskWeights[tasksLeft - tasksToShare]), tasksToShare, MPI_INT, rankOfClient, TASKS_SENDING_TAG, MPI_COMM_WORLD);
            tasksLeft -= tasksToShare;
            pthread_mutex_unlock(&mutex);
        }
        else{
            pthread_mutex_unlock(&mutex);
            tasksToShare = 0;
            MPI_Send(&tasksToShare, 1, MPI_INT, rankOfClient, TASK_NUMBER_REPLY_TAG, MPI_COMM_WORLD);
        }
    }
}


int main(int argc, char* argv[]) {

    int threadSupportLevel;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &threadSupportLevel);
    if(threadSupportLevel != MPI_THREAD_MULTIPLE){
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rankOfProc);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcs);

    pthread_mutex_init(&mutex, NULL);
    pthread_attr_t threadAttributes;
    pthread_attr_init(&threadAttributes);
    pthread_attr_setdetachstate(&threadAttributes, PTHREAD_CREATE_JOINABLE);

    pthread_create(&threads[0], &threadAttributes, managersWork, NULL);
    pthread_create(&threads[1], &threadAttributes, executorsWork, NULL);

    pthread_join(threads[0], NULL);
    pthread_join(threads[1], NULL);
    pthread_attr_destroy(&threadAttributes);
    pthread_mutex_destroy(&mutex);

    MPI_Finalize();
    return EXIT_SUCCESS;
}