

/* strtok example by mind@metalshell.com
 *
 * This is an example on string tokenizing
 *
 * 02/19/2002
 *
 * http://www.metalshell.com
 *
 */
 
#include <stdio.h>
#include <string.h>
 
/* Get last word divided by "/"  CB  3-MAY-2010  */
char* token(char *str)
{
        int x = 1;
        char *str1;
        char *tok;
 
        /* print what we have so far  */
        /* printf("String: %s\n", str);   */
 
        /* extract first string from string sequence */
        str1 = strtok(str, "/");
 
        /* print first string after tokenized */
        /* printf("%i: %s\n", x, str1);  */
        tok = str1;
 
        /* loop until finished */
        while (str1 != NULL)
        {
                /* extract string from string sequence */
                str1 = strtok(NULL, "/");
 
                /* check if there is nothing else to extract */
                if (str1 == NULL) break;
               
                tok = str1;
                /* print string after tokenized  */
                /* printf("%i: %s\n", x, str1);  */
                x++;
        }
        /* printf("%i: %s\n", x, tok);   */
 
        return tok;
 
}


