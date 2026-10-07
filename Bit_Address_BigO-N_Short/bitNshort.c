#include <stdio.h>
#include <stdlib.h>

struct Bitdata
{
    int data;
    struct Bitdata* next;

    /* data */
};

struct locData{
    int data;
    int Max;
    struct locData* next;
};

int getBit(int num, int bitPosition) {
    return (num >> bitPosition) & 1; // Shift right and mask with 1 to get the bit
}

int isNumInArr(int num, int arr[], int size) {
    for(int i = 0; i < size; i++) {
        if(arr[i] == num) {
            return 1; // Number found in the array
        }
    }
    return 0; // Number not found in the array
}

void changeNumLocinArr(int arr[],int size,int n,int newLoc){
    int arrtool[size+1];
    int cnt=0;
    for(int i=0; i<size; i++){
        if(cnt==0 && i == newLoc){
            arrtool[newLoc] = arr[n];
            cnt = 1;
        }
        else if(cnt==1){
            arrtool[i]=arr[i-1];
        }
        else if(cnt==1 && i==n) {
            arrtool[i] = arr[i-1];
            cnt=2;
        }
        else arrtool[i] = arr[i];
    }
}

int main(){
 
    int n=2;
    
    int arr[n];

    for(int i = 0; i < n; i++){
        scanf("%d", &arr[i]);
    }


    struct Bitdata bitArr[32];

    for(int i=0; i<32; i++){
        bitArr[i].data = 0; // Initialize the bit array
        bitArr[i].next = NULL; // Initialize the next pointer
        struct Bitdata* BitdataTool = &bitArr[i]; // Create a pointer to the current bit array element
        for(int j=0; j<n; j++){
            int tool = getBit(arr[j], i);
            if(tool == 1){
                BitdataTool->next = (struct Bitdata*)malloc(sizeof(struct Bitdata)); // Allocate memory for the next node
                BitdataTool->next->data = j+1; // Store the original number in the next node
                BitdataTool->next->next = NULL; // Initialize the next pointer of the new node
                BitdataTool = BitdataTool->next; // Move the pointer to the new node
            }
        }
    }
    
    for(int i=0; i<32; i++){
        struct Bitdata* tool = bitArr[i].next;
        printf("bit %d: ", i);
        while(tool != NULL){
            printf("%d ", tool->data);
            tool = tool->next;
        }
        printf("\n");
    }

    int arrshort[n];
    struct locData shortdata[32];
    struct locData* locOfShortData[n];

    int index = 0;
    for(int i=31; i>=0; i--){
        int indextool = index;

        struct Bitdata* tool = bitArr[i].next;

        struct locData* shortTool = &shortdata[i];
        shortTool->data = indextool; // Initialize the shortdata array
        shortTool->Max = indextool; // Initialize the Max value in the shortdata array
        shortTool->next = NULL; // Initialize the next pointer in the shortdata array

        int indextool2=indextool;
        int indextool2Num=0;

        while(tool != NULL){
            if(!isNumInArr(tool->data, arrshort, index)){ // Check if the number is already in the short array
                arrshort[index] = tool->data; // Store the original number in the short array
                shortdata[i].next = (struct locData*)malloc(sizeof(struct locData)); // Allocate memory for the next node in the shortdata array
                
                shortTool = shortTool->next; // Move to the next node in the shortdata array
                shortTool->data = tool->data; // Store the original number in the shortdata array
                shortTool->Max = indextool; // Store the bit position in the shortdata array
                shortTool->next = NULL; // Initialize the next pointer
                locOfShortData[tool->data-1]= shortTool; // Store the location of the shortdata node in the locOfShortData array
                
                index++;
            }
            else{
                changeNumLocinArr(arrshort, index, tool->data-1, indextool2); // Change the location of the number in the short array
                indextool2++;
                indextool2Num++;
            }
            tool = tool->next; // Move to the next node
        }

        if(indextool2Num>0){
            for(int j=indextool2; j<index;j++){
                locOfShortData[arrshort[j]-1]->Max+=indextool2Num; // Update the Max value in the shortdata array
            }
        }
    }

    for(int i=0; i<n; i++){
        printf("%d ", arrshort[i]);
    }
    return 0;
}