/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 **/

#include <stdlib.h>

#include "plasma_context.h"
#include "plasma_internal.h"

static int max_contexts = 1024;
static int num_contexts = 0;

plasma_context_map_t *context_map = NULL;
pthread_mutex_t context_map_lock = PTHREAD_MUTEX_INITIALIZER;

/***************************************************************************//**
    @ingroup plasma_init
    Initializes PLASMA, allocating its context.
*/
int PLASMA_Init()
{
    pthread_mutex_lock(&context_map_lock);

    // Allocate context map if NULL.
    if (context_map == NULL) {
        context_map =
            (plasma_context_map_t*)calloc(max_contexts,
                                          sizeof(plasma_context_map_t));
        if (context_map == NULL) {
            pthread_mutex_unlock(&context_map_lock);
            plasma_error("calloc() failed");
            return PlasmaErrorOutOfMemory;
        }
    }
    pthread_mutex_unlock(&context_map_lock);

    plasma_context_attach();
    return PlasmaSuccess;
}

/***************************************************************************//**
    @ingroup plasma_init
    Finalizes PLASMA, freeing its context.
*/
int PLASMA_Finalize()
{
    plasma_context_detach();
    return PlasmaSuccess;
}

/******************************************************************************/
int PLASMA_Set(plasma_enum_t param, int value)
{
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_error("PLASMA not initialized");
        return PlasmaErrorNotInitialized;
    }
    switch (param) {
    case PLASMA_TILE_SIZE:
        if (value <= 0) {
            plasma_error("invalid tile size");
            return PlasmaErrorIllegalValue;
        }
        plasma->nb = value;
        break;
    case PLASMA_INNER_BLOCK_SIZE:
        if (value <= 0) {
            plasma_error("invalid inner block size");
            return PlasmaErrorIllegalValue;
        }
        plasma->ib = value;
        break;
    default:
        plasma_error("Unknown parameter");
        return PlasmaErrorIllegalValue;
    }
    return PlasmaSuccess;
}

/******************************************************************************/
int PLASMA_Get(plasma_enum_t param, int *value)
{
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA not initialized");
        return PlasmaErrorNotInitialized;
    }
    switch (param) {
    case PLASMA_TILE_SIZE:
        *value = plasma->nb;
        return PlasmaSuccess;
        break;
    case PLASMA_INNER_BLOCK_SIZE:
        *value = plasma->ib;
        return PlasmaSuccess;
        break;
    default:
        plasma_error("Unknown parameter");
        return PlasmaErrorIllegalValue;
    }
    return PlasmaSuccess;
}

/******************************************************************************/
int plasma_context_attach()
{
    pthread_mutex_lock(&context_map_lock);

    // Reallocate context map if out of space.
    if (num_contexts == max_contexts-1) {
        max_contexts *= 2;
        context_map = (plasma_context_map_t*) realloc(
            &context_map, max_contexts*sizeof(plasma_context_map_t));
        if (context_map == NULL) {
            pthread_mutex_unlock(&context_map_lock);
            plasma_error("realloc() failed");
            return PlasmaErrorOutOfMemory;
        }
    }
    // Create the context.
    plasma_context_t *context;
    context = (plasma_context_t*)malloc(sizeof(plasma_context_t));
    if (context == NULL) {
        pthread_mutex_unlock(&context_map_lock);
        plasma_error("malloc() failed");
        return PlasmaErrorOutOfMemory;
    }
    // Initialize the context.
    plasma_context_init(context);

    // Find and empty slot and insert the context.
    for (int i = 0; i < max_contexts; i++) {
        if (context_map[i].context == NULL) {
            context_map[i].context = context;
            context_map[i].thread_id = pthread_self();
            num_contexts++;
            pthread_mutex_unlock(&context_map_lock);
            return PlasmaSuccess;
        }
    }
    // This should never happen.
    pthread_mutex_unlock(&context_map_lock);
    plasma_error("empty slot not found");
    return PlasmaErrorInternal;
}

/******************************************************************************/
int plasma_context_detach()
{
    pthread_mutex_lock(&context_map_lock);

    // Find the thread and remove its context.
    for (int i = 0; i < max_contexts; i++) {
        if (context_map[i].context != NULL &&
            pthread_equal(context_map[i].thread_id, pthread_self())) {

            free(context_map[i].context);
            context_map[i].context = NULL;
            num_contexts--;
            pthread_mutex_unlock(&context_map_lock);
            return PlasmaSuccess;
        }
    }
    pthread_mutex_unlock(&context_map_lock);
    plasma_error("context not found");
    return PlasmaErrorInternal;
}

/******************************************************************************/
plasma_context_t *plasma_context_self()
{
    pthread_mutex_lock(&context_map_lock);

    // Find the thread and return its context.
    for (int i = 0; i < max_contexts; i++) {
        if (context_map[i].context != NULL &&
            pthread_equal(context_map[i].thread_id, pthread_self())) {

            pthread_mutex_unlock(&context_map_lock);
            return context_map[i].context;
        }
    }
    pthread_mutex_unlock(&context_map_lock);
    plasma_error("context not found");
    return NULL;
}

/******************************************************************************/
void plasma_context_init(plasma_context_t *context)
{
    context->nb = 256;
    context->ib = 64;
    context->translation = PLASMA_OUTOFPLACE;
}
