#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <regex.h>

/* -----------------------------------------------
 * Searches through a directory, and returns
 * frequencies in the pN file names (sorted).
 *
 * Started 5.25.2016, Austin Sousa -- asousa@stanf
 * ----------------------------------------------- 
 */





typedef struct node {
    int freq;
    float L;
    struct node *ptr;
} node;


// Prototypes
node* insert(node* head, int freq, float L);
void free_list(node *head);
void print_list(node* head);
int * get_freqs_at(char *directory, int *num_freqs);

// Search a directory for pN, pS files
// Return a sorted list of frequencies to work with. BAM.


// int main(void)

// {
//     int *freqs;
//     int num_freqs;
//     int i;

//     freqs = get_freqs_at("./out_matching_mycrossings", &num_freqs);

//     printf("Recipt: %d\n",num_freqs);

//     for (i=0; i<num_freqs; i++) {
//         printf("i: %d, f: %d\n", i,freqs[i]);
//     }

// }




/* -----------------------------------------------
 * Returns a pointer to an array of integers.
 * Assigns the total length of the array to num_freqs.
 * -----------------------------------------------
 */
int * get_freqs_at(char *directory, int *num_freqs)
{
  DIR *dp;
  struct dirent *ep;     
  dp = opendir (directory);
  int n=0;
  int s=0;
  int i=0;
  float L_N, L_S;
  int freq_N, freq_S;
  // char *strname;
  int *freq_output;


  node *list_N = NULL;
  node *list_S = NULL;
  node *p;

  if (dp != NULL)
  {
    while (ep = readdir (dp)) {
      if (sscanf(ep->d_name, "pN%d_%g.dat", &freq_N, &L_N)) {
        // printf("N: %d, %g\n",freq_N, L_N);
        list_N = insert(list_N, freq_N, L_N);
        n++;
      }

      if (sscanf(ep->d_name, "pS%d_%g.dat", &freq_S, &L_S)) {
        // printf("S: %d, %g\n",freq_S, L_S);
        list_S = insert(list_S, freq_S, L_S);
        s++;
      }
    }

    closedir (dp);

    printf("n is: %d, s is: %d\n",n,s);

    // print_list(list_N);
    // print_list(list_S);

    freq_output = (int*)calloc(n,sizeof(int));

    p = list_N;
    for (i=0; i<n; i++) {
        // while(p) {
            // printf("%d ", p->freq);
            freq_output[i] = p->freq;
            // printf("%d\n",freq_output[i]);
            p = p->ptr;
        // }
    }
  
    *num_freqs = n;

    free_list(list_N);
    free_list(list_S);

    return freq_output;    
    // for (i=0; i<n; i++) {
    //     printf("%d\n",freq_output[i]);
    // }


  //   for (i=0; i<n; i++) {
  //       printf("i: %d, f: %d\n",i,freq_output[i]);
  //   }

  // printf("Got %d frequencies\n",n);
  }


  else
    perror ("Couldn't open the directory");

  // return 0;
}



// -----------------------------------------------
// A simple linked-list.
// -----------------------------------------------
node* insert(node* head, int freq, float L) {
    node *temp, *prev, *next;


    temp = (node*)malloc(sizeof(node));
    temp->freq = freq;
    temp->L = L;
    temp->ptr = NULL;

    if (!head) {
        head=temp;
    } else {
        prev = NULL;
        next = head;
        while(next && next->freq <= freq){
            prev = next;
            next = next->ptr;
        }
        if(!next){
            prev->ptr = temp;
        } else {
            if(prev) {
                temp->ptr = prev->ptr;
                prev-> ptr = temp;
            } else {
                temp->ptr = head;
                head = temp;
            }            
        }   
    }
    return head;
}

// Print the list contents
void print_list(node* head) {
    node * p;

    p = head;
    while(p) {
        printf("%d ", p->freq);
        p = p->ptr;
    }
}

// Delete the list, free up memory.
void free_list(node *head) {
    node *prev = head;
    node *cur = head;
    while(cur) {
        prev = cur;
        cur = prev->ptr;
        free(prev);
    }       
}