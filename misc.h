#pragma once

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/*
 * Return if a path is a exsiting directory.
 */
inline
bool is_dir(char *path) {
    struct stat sb = {};
    stat(path, &sb);
    return (sb.st_mode & S_IFMT) == S_IFDIR;
}

/*
 * Return if a path is a exsiting file.
 */
inline
bool is_file(char *path) {
    struct stat sb = {};
    stat(path, &sb);
    return (sb.st_mode & S_IFMT) == S_IFREG;
}
