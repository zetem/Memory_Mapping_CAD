#pragma once
#ifndef HEAD_H
#define HEAD_H

#include <string>
using namespace std;
//LUTRAM : 64x10bit or 32x20bit
#define N                   10 //# of LUTs in each block

#define LB_COUNT_LUTRAM     2 // 50 percent so for every 2 blocks one can be LUTRAM
#define MAX_SIZE_LUTRAM     640 //1024*128
#define MAX_WIDTH_LUTRAM    20  //for depth 64 it is half

#define LB_COUNT_8K         10
#define MAX_WIDTH_8K        32  //for tru dual it is half
#define MAX_SIZE_8K         8192 //1024*8

#define LB_COUNT_128K       300
#define MAX_WIDTH_128K      128 //for tru dual it is half
#define MAX_SIZE_128K       131072 //1024*128

#define MAX_SERIAL          16

#define RAM_TYPES_COUNT     4 
enum    RAMtype             {ROM , SinglePort, SimpleDualPort, TrueDualPort};
#define MAX_PHRAM_TYPES_COUNT   3

#define LB_AREA_LUT         40000
#define LB_AREA             35000

//enum    PHRAMtype           {LUTRAM, M8K, M128K};
enum    PHRAMtype           {TYPE0, TYPE1, TYPE2};
enum    DIR                 {MOVEBIGGER, MOVESMALLER};

int PHRAM_TYPES_COUNT = 0;
bool STRATIX = false;
bool ENABLE_COUT = false;

class RAM;
class Circuit;
RAMtype RAMtype_num(string s);

class Bad_Data{}; 

class PhRAM{
public: 
    //PhRAM();
    PhRAM(int t);
    PhRAM(int s, int w , int r, int t): max_size(s), max_width(w), ratio(r), type(t){area = (9000 + 5*s + 90*sqrt(s) + 600*2*w);}
    virtual int cost_to_map(RAM* r, int &parallel, int &serial, int& width) const;
    virtual int cost_to_map(int d , int w , bool td, int &parallel, int &serial, int& width) const;
    //virtual void unset_LUTRAM() = 0;
    //PhRAM* create_phRAM();
    int type;
    int max_width;
    int max_size;
    int ratio;
    long area;
};

class LutRAM : public PhRAM {
public: 
    LutRAM(): PhRAM(0) {}
    //void unset_LUTRAM() {support_LUTRAM = false;} 
    int cost_to_map(RAM* r, int &parallel, int &serial, int& width) const;
    int cost_to_map(int d , int w , bool td, int &parallel, int &serial, int &width) const;
};


/////////////////////////////////////////////////////////////////
PhRAM* phRAMs[MAX_PHRAM_TYPES_COUNT];
//vector<PhRAM*> phRAMs;
//=============================================================================//

class Physical_RAM_Bank{
public:
    Physical_RAM_Bank( int i, int w, int s, int p, int r, Circuit* c, int t, int add): id(i), width(w), serial(s), parallel(p), type(t), logic_added(add), mapped_logical_index(r), top_circuit(c){}
    Physical_RAM_Bank* copy();
    void modify(int s, int p, int w, int t,int lw, bool td);
    void modify(Physical_RAM_Bank* b);
protected:
    int id; //should be the same as logical ID cause using a bank of RAMs for each logical RAM
    int width;
    int serial; //serial
    int parallel; //parallel
    int logic_added;
    int mapped_logical_index; //logical RAM
    //int type;
    int type;
    Circuit *top_circuit;
    friend ostream& operator<<(ostream&, const Physical_RAM_Bank*); 
    friend class Circuit;
};

class RAM {
public:
    RAM(int i, string t, int d, int w , Circuit* c)
        : id(i), type(RAMtype_num(t)) , depth(d), width(w), mapped(false), top_circuit(c) {}

    //virtual void play_turn();
    bool is_mapped() {return mapped;}
    void map_RAM() {mapped = true;}
    bool less(RAM* rc, int b);
    bool equal(RAM* rc); 
    void map_lowest_individual_cost();   
    int RAM_size();

protected:
    int id; 
    RAMtype type; // r for Rom , s for single port , d for simple dual port, t for true dual port 
    int depth;
    int width;
    int get_id(){return id;} // logical ID
    bool mapped; 
    Physical_RAM_Bank* mapped_RAM;
    Circuit *top_circuit; //! should remove no use 

    friend ostream& operator<<(ostream&, const RAM&);
    friend ostream& operator<<(ostream&, const RAM*);
    friend ostream& operator<<(ostream&, const Physical_RAM_Bank* );
    friend class Circuit;
    friend class PhRAM;
    friend class LutRAM;
    //friend class Physical_RAM_Bank;
};


//=======================================================================//
class Circuit {
public:
    Circuit(int id, int lbc):circuit_id(id), logic_block_count(lbc){}
    ~Circuit();
    void add_logical_RAM(RAM* new_RAM);
    void update_available_resources();
    void update_required_resources();
    int  update_max_resources(int &index);
    void optmize_area();
    void update_area();
    long get_geo_area(){return geo_area;}
    bool modify_mapping(int t, int small_big);
    int get_id() {return circuit_id;}
    vector<RAM*> logical_RAMs; // made it public for quick sort to change
    vector<Physical_RAM_Bank*> physical_RAMs;
    vector<int> sorted_logical_RAMs[MAX_PHRAM_TYPES_COUNT][2]; //2 for one MOVEBIGGER and one MOVESMALLER 
private:
    int circuit_id;
    int logic_block_count;
    int available_resources[MAX_PHRAM_TYPES_COUNT];
    int max_available_resources[MAX_PHRAM_TYPES_COUNT];
    int required_resources[MAX_PHRAM_TYPES_COUNT+1]; //+1 for added LUT by dec and mux of memories
    int area; //total logic block count
    long geo_area;
    int limiting_rsc;
    vector<int> unmapped_tdual;  // if tru dual push back the id, true dual cannot use LUTRAM

    friend ostream& operator<<(ostream&, const Circuit&);
    friend ostream& operator<<(ostream& out, const Circuit*); 
};


#endif
