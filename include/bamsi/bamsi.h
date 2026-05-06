#ifndef BAMSI_BAMSI_H
#define BAMSI_BAMSI_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct bamsi_index bamsi_index_t;

typedef enum bamsi_status_t {
    BAMSI_STATUS_OK = 0,
    BAMSI_STATUS_INVALID_ARGUMENT = 1,
    BAMSI_STATUS_NOT_IMPLEMENTED_V1 = 2
} bamsi_status_t;

#ifdef __cplusplus
}
#endif

#endif  // BAMSI_BAMSI_H
