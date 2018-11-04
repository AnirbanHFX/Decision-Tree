#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ROWS 24
#define TESTROWS 4
#define COLS 8

typedef struct _node {   //Node for tree data structure
   int attribute; 
   /*
   *attribute = -2 indicates FALSE
   *attribute = -1 indicates TRUE
   *attribute >= 0 indicates attribute based on which further nodes have been split
   */
   struct _node *left;
   struct _node *right;
} node;

typedef struct _table {   //Structure to store truth table
   int avail;   //Checks whether row is available. Used to restrict truthtable to a subset
   int *x;
   int result;  //Result stores value of y
} tab;

double gain(tab *table, int index)   //Information gain function

{

   int i;
   double total = 0, t=0, f=0, t0=0, t1=0, f0=0, f1=0;

   /*
   *total - total available rows
   *t - number of available rows with y=1
   *f - number of available rows with y=0
   *t0 - number of available rows with y=1 and x[index] = 0
   *t1 - number of available rows with y=1 and x[index] = 1
   *f0 - number of available rows with y=0 and x[index] = 0
   *f1 - number of available rows with y=0 and x[index] = 1
   */

   double ei, e0, e1;

   /*
   *ei - Initial entropy before splitting
   *e0 - Entropy of all available rows with y=0
   *e1 - Entropy of all available rows with y=1
   */

   for (i=0; i<ROWS; i++ ) {

      if(table[i].avail == 1) {

         total += 1;

         if(table[i].result == 1) {
            t+=1;
            if(table[i].x[index] == 1) t1 += 1;
            else t0 += 1;
         }
         else if(table[i].result == 0) {
            f+=1;
            if(table[i].x[index] == 1) f1 += 1;
            else f0 += 1;
         }

      }

   }

   if (total == 0) return 0;

   if (t==0 || f==0)
      ei = 0;
   else
      ei = -(t/total)*log(t/total) - (f/total)*log(f/total);

   if (t0==0 || f0==0)
      e0 = 0;
   else
      e0 = -(t0/(t0+f0))*log(t0/(t0+f0)) - (f0/(t0+f0))*log(f0/(t0+f0));

   if (t1==0 || f1==0)
      e1 = 0;
   else
      e1 = -(t1/(t1+f1))*log(t1/(t1+f1)) - (f1/(t1+f1))*log(f1/(t1+f1));

   return(ei - ((t0+f0)/total)*e0 - ((t1+f1)/total)*e1)/log(2.0);  //Returns information gain

}





node* createTree(tab *table, int *att) {   //Creates decision tree

   node *leaf = (node *) malloc(sizeof(node));
   int index;   //Stores index of attribute which furnishes maximum information gain on splitting

   int attl[COLS], attr[COLS];

   /*
   *att[i] = 1 if ith attribute has not been used, 0 otherwise
   *attl is passed to the left child node after required changes
   *attr is passed to the right child node after required changes
   */

   for(int i=0; i<COLS; i++) {  //Initializes attl and attr with att
      attl[i] = att[i];
      attr[i] = att[i];
   }

   int total=0, t=0, f=0;

   /*
   *total - total available rows
   *t - number of available rows with y=1
   *f - number of available rows with y=0
   */

   for(int i=0; i<ROWS; i++) {   //Counts total, t, f

      if(table[i].avail == 1) {
         total++;
         if(table[i].result == 1) t++;
         else if(table[i].result == 0) f++;
      }

   }

   if (t==total) {
      leaf->attribute = -1;
      return leaf;   //No further splitting required
   }

   else if (f==total) {
      leaf->attribute = -2;
      return leaf;   //No further splitting required
   }

   int over = 1;   //Over is initialized to 1 indicating no attributes are left to be split on

   for(int i=0; i<COLS; i++) {
      if(att[i] != 0) {    //Checks whether any attribute is left
         over = 0;   //If an unused attribute is found, over is set to 0
         break;
      }
   }

   if (over == 1) {   //If no attributes are left, further splitting cannot be done
      if(t>=f) leaf->attribute = -1;
      else leaf->attribute = -2;
      return leaf;
   }

   double highgain = -1;   //Highest information gain initialized to negative value

   for (int i=0; i<COLS; i++) {

      if(att[i]!=0&&(gain(table, i)) > highgain) {
         highgain = gain(table, i);
         index = i;   //Stores index of attribute furnishing highest information gain
      }

   }

   tab *tableleft = (tab *) malloc(ROWS*sizeof(tab));  //x[index] is 0
   tab *tableright = (tab *) malloc(ROWS*sizeof(tab));  //x[index] is 1

   for(int i=0; i<ROWS; i++) {   //Creates two truthtables that are subsets of the original, to be passed on to child nodes
      tableleft[i] = table[i];
      tableright[i] = table[i];
      if(table[i].x[index] == 1) tableleft[i].avail = 0;
      else if (table[i].x[index] == 0) tableright[i].avail = 0;
   }

   int leftc = 0, rightc = 0;

   /*
   *leftc - Number of available rows in tableleft
   *rightc - Number of available rows in tableright
   */

   for(int i=0; i<ROWS; i++) {
      if(tableleft[i].avail == 1) leftc++;
      if(tableright[i].avail == 1) rightc++;
   }

   leaf->attribute = index;   //Assigns attribute on which data is to be split

   leaf->left = (node *) malloc(sizeof(node));
   leaf->right = (node *) malloc(sizeof(node));

   if (leftc == 0) {   //tableleft is empty, so assigns most common y from parent node to left child
      if (t>=f)
         leaf->left->attribute = -1;
      else
         leaf->left->attribute = -2; 
   }

   else {
      attl[index] = 0;   //Removes indexed attribute from attl
      leaf->left = createTree(tableleft, attl);   //Recursion on left child
   }

   if (rightc == 0) {   //tableright is empty, so assigns most common y from parent node to right child
      if (t>=f)
         leaf->right->attribute = -1;
      else
         leaf->right->attribute = -2; 
   }

   else {
      attr[index] = 0;   //Removes indexed attribute from attr
      leaf->right = createTree(tableright, attr);   //Recursion on right child
   }

   return leaf;

}

int findTree(node *bud, int *instance)   //Recursively classifies testcases

{

   if(bud->attribute < 0) return (bud->attribute + 2);

   else if(instance[bud->attribute] == 0) bud = bud->left;
   
   else if(instance[bud->attribute] == 1) bud = bud->right;

   return findTree(bud, instance);

}

int main()

{

   int att[COLS];
   FILE *fp, *fp2, *fpw;
   int i, j;
   tab *table = (tab *) malloc(ROWS*sizeof(tab));

   int **testdata = (int **) malloc(TESTROWS*sizeof(int *));

   fp = fopen("data2.csv", "r");

   if (fp == NULL) {
      printf("data2.csv not found\n");
      return 1;
   }

   for(j=0; j<COLS; j++) att[j] = 1;

   for (i=0; i<ROWS; i++) {

      table[i].avail = 1;
      table[i].x = (int *) malloc(COLS*sizeof(int));

      for (j=0; j<COLS; j++) {

         fscanf(fp, " %d,", &table[i].x[j]);

      }

      fscanf(fp, " %d,", &table[i].result);

   }

   node *head = (node *) malloc(sizeof(node));

   head = createTree(table, att);  //Creates tree

   fclose(fp);

   fp2 = fopen("test2.csv", "r");

   if (fp2 == NULL) {
      printf("test.csv not found\n");
      return 1;
   }

   for(i=0; i<TESTROWS; i++) {

      testdata[i] = (int *) malloc(COLS*sizeof(int));

      for(j=0; j<COLS; j++) {

         fscanf(fp2, " %d,", &testdata[i][j]);

      }

   }

   fclose(fp2);

   fpw = fopen("Decision_tree.out", "w");

   if(fpw == NULL) {
      printf("Unable to create output file\n");

      for(i=0; i<TESTROWS; i++) 
         printf("%d ", findTree(head, testdata[i]));

      printf("\n");

      return 1;
   }

   for(i=0; i<TESTROWS; i++) {
      printf("%d ", findTree(head, testdata[i]));
      fprintf(fpw, "%d ", findTree(head, testdata[i]));
   }

   fclose(fpw);

   printf("\n");   

   return 0;

}
