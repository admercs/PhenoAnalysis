/* Copyright 2008-2013 Colorado State University and Blue Marble Research
   Authors: Matthew Bishop (Colorado State University) and Reto Stockli (Blue Marble Research) 
   
   This file is part of <INSERT_PACKAGE_NAME> and was created within the
   EUMETSAT Satellite Application Facility on Climate Monitoring (CM SAF)
   project by MeteoSwiss.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>. */

/* NOT THREAD SAFE. */

#include <assert.h>
#include <glob.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_FILEPAT_LEN 512

static glob_t matches;
static int current_file = -1;

/* Matches files based on the given file pattern.  Returns the number of
 * files that were found.
 */
int match_file_init(int *len, char *filepat)
{
    char cpat[MAX_FILEPAT_LEN];
    int i, ret;

    /* Get filepat into a proper null-terminated C string. */
    assert(*len < MAX_FILEPAT_LEN);
    for (i = 0; i < *len; i++)
        cpat[i] = filepat[i];
    cpat[i] = '\0';

    ret = glob(cpat, 0, NULL, &matches);
    switch (ret) {
    case 0: /* no errors */
        break;
    case GLOB_NOMATCH: /* no matches found */
        return 0;
    case GLOB_NOSPACE: /* ran out of memory */
        fprintf(stderr, "glob(): ran out of memory\n");
        abort();
    case GLOB_ABORTED: /* read error */
    default:
        fprintf(stderr, "glob(): unexpected error %i\n", ret);
        abort();
    }
    current_file = 0;
    return (int)matches.gl_pathc;
}

/* Returns the next matched file from the call to fl_match_files().  Only
 * call as many times as the return value from fl_match_files() indicates.
 * file: out parameter containing the next file's name.
 * len:  implicit parameter containing the length of file.
 */
void match_file_next(int *len, char *file)
{
    char *s = matches.gl_pathv[current_file++];
    int i = 0;

    while (i < *len && s[i] != '\0') {
        file[i] = s[i];
        i++;
    }
    while (i < *len)
        file[i++] = ' ';

}

/* Frees the resources used from the call to fl_match_files().  Call this
 * once after every call to fl_match_files() but only after you are done
 * calling fl_next_file().
 */
void match_file_exit(void)
{
    globfree(&matches);
    current_file = -1;
}


/* simple test code */
/*
int main(int argc, char **argv)
{
    char filebuf[1024];
    int i, n;

    n = match_file_init(argv[argc - 1], (int *) strlen(argv[argc - 1]));
    for (i = 0; i < n; i++) {
        match_file_next(filebuf, (int *) sizeof(filebuf) - 1);
        printf("%i: %s\n", i, filebuf);
    }
    match_file_exit();

    return 0;
} 
*/
