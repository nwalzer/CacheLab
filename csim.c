#include "cachelab.h"
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <getopt.h>

/*
	Nathan Walzer - nwalzer
	Lucas Varella - lnvarella
	Group: lnvarella-nwalzer
*/
#define VALID 'T'
#define INVALID 'F'

struct block {
    char valid;
    unsigned long tag;
    int LRU;
    //for the purposes of this assignment we can ignore the bytes that would be stored
};

struct block** alloCache(unsigned int s, int e){//allocates the cache to the correct size
    struct block** tempCache = (struct block**) calloc(1<<s, sizeof(struct block));
    if(tempCache == NULL) return NULL;
    for(int i = 0; i < 1<<s; i++){
	tempCache[i] = (struct block*) calloc(e, sizeof(struct block));
        if(tempCache[i] == NULL) return NULL;
	for(int j = 0; j < e; j++){
	    tempCache[i][j].valid = INVALID;
	    tempCache[i][j].LRU = 0;
	}
    }
    return tempCache;
}

//if the given set contains a valid line with the given tag return "true" (1)
int isHit(struct block** cache, unsigned int set, unsigned long tag, int lines){
    for(int i = 0; i < lines; i++){
	if(cache[set][i].valid == VALID && cache[set][i].tag == tag){
	    cache[set][i].LRU = 0;
	    return 1;
	}
    }
    return 0;
}

//return the index of the first invalid line in a set, otherwise return -1
int anyInvalid(struct block** cache, unsigned int set, int lines){
    for(int i = 0; i < lines; i++){
	if(cache[set][i].valid == INVALID) return i;
    }
    return -1;
}

//place the given block at the given location, overwriting the previous data
void place(struct block** cache, unsigned int set, int idx, unsigned long tag){
    cache[set][idx].valid = VALID;
    cache[set][idx].tag = tag;
    cache[set][idx].LRU = 0;
    return;
}

//evict the least recently used block and set its LRU count to 0
void evict(struct block** cache, unsigned int set, int lines, unsigned long tag){
    int maxLRU = 0;
    int mLRUIdx = 0;
    for(int i = 0; i < lines; i++){
	if(cache[set][i].LRU > maxLRU){
	    maxLRU = cache[set][i].LRU;
	    mLRUIdx = i;
	}
    }
    place(cache, set, mLRUIdx, tag);
}

//go through a given set and increment all of the LRUs
void incLRU(struct block** cache, int set, int lines){
    for(int i = 0; i < lines; i++) cache[set][i].LRU++;
}

int main(int argc, char** argv){
    int s;
    int e;
    int b;
    int opt;
    int hits = 0;
    int miss = 0;
    int evic = 0;
    unsigned int set;
    struct block** cache;
    FILE *t;
    
    while((opt = getopt(argc, argv, "s:E:b:t:")) != -1){
	switch(opt){
	case 's':
	    s = atoi(optarg);
	    break;
	case 'E':
	    e = atoi(optarg);
	    break;
	case 'b':
	    b = atoi(optarg);
	    break;
	case 't':
	    t = fopen(optarg, "r");
		printf("%s", optarg);
	    break;
	case '?': //if we get unexpected input abort program
	    return 0;
	default:
	    break;
	}
    }
    //s = atoi(argv[2]);
    //e = atoi(argv[4]);
    //b = atoi(argv[6]);
    //t = fopen(argv[8], "r");
    if(t == NULL) return 0; //if the file didn't open exit the program
    cache = alloCache(s, e);
    if(cache == NULL) return 0; //if the cache wasn't allocated exit the program

    char* op = malloc(8);
    unsigned long addr;
    char* after = malloc(8);
    int invalIdx;
    while (fscanf(t, "%s", op) != EOF){ //get the type of instruction (I, M, S, L)
        fscanf(t, "%lx", &addr); //get the address
    	fscanf(t, "%s", after); //unused space after the address
	if(op[0] == 'I') continue; //if it's an instruction argument then ignore
	addr = addr>>b;//ignore the offset bits
	set = addr & ~(0x7FFFFFFFFFFFFFFFL<<s);//isolate the set bits
	addr = addr>>s;//addr now becomes tag bits
	incLRU(cache, set, e);//increment all LRUs
	if(op[0] == 'M') hits++;//"M" always guarentees at least one hit
	if(isHit(cache, set, addr, e)){//if it is a hit then increment hits and restart the loop
	    hits++;
	    continue;
	}
	miss++;//if not a hit, then inc miss
	invalIdx = anyInvalid(cache, set, e);//Place block in invalid line before evicting other lines
	if(invalIdx != -1){
	    //if invalIdx returns an index, place this block in that postiion
	    place(cache, set, invalIdx, addr);
	} else {
	    //otherwise we must evict the highest LRU
	    evic++;
	    evict(cache, set, e, addr);
	}
    }

    printSummary(hits, miss, evic);
    return 0;
}
