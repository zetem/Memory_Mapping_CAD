#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <stdlib.h>     /* atoi */
#include <ctime>
#include <fstream>
#include <sstream>
#include <math.h>
#include<time.h>
//#include <algorithm>
//#include <limits>
#include "HEAD.h"
using namespace std;


RAMtype RAMtype_num(string s)
{   
    if (s=="ROM")
        return ROM;
    else if (s=="SinglePort")
        return SinglePort;
    else if (s=="SimpleDualPort")
        return SimpleDualPort;
    else if (s=="TrueDualPort")
        return TrueDualPort;
}

string RAMtype_name(RAMtype t)
{
    string RAM_name[RAM_TYPES_COUNT] = {"ROM", "SinglePort", "SimpleDualPort" ,"TrueDualPort"};
    return RAM_name[t];
}
string PHRAMtype_name(int t)
{
    string RAM_name_d[MAX_PHRAM_TYPES_COUNT] = {"LUTRAM", "8K", "128K"};
    string RAM_name[MAX_PHRAM_TYPES_COUNT]   = {"TYPE0", "TYPE1", "TYPE2"};
    if (STRATIX)
        return RAM_name_d[t];
    return RAM_name[t];
}


PhRAM::PhRAM(int t){
    type = t; 
    switch (t){
        case 0://LUTRAM:
            max_size = MAX_SIZE_LUTRAM;
            max_width = MAX_WIDTH_LUTRAM;
            ratio = LB_COUNT_LUTRAM;
            break;
        case 1://M8K:
            max_size = MAX_SIZE_8K;
            max_width = MAX_WIDTH_8K;
            ratio = LB_COUNT_8K;
            break;
        case 2://M128K:
            max_size = MAX_SIZE_128K;
            max_width = MAX_WIDTH_128K;
            ratio = LB_COUNT_128K;
            break;
    }
    area = (9000 + 5*max_size + 90*sqrt(max_size) + 600*2*max_width);
}
//////////////////////////////////////////////////////////////////
int added_logic_count(int s,int lw, bool td)
{       
    int mux_lut = 0;
    int dec_lut = 0; // if numner of inputs is less than 5( 2 5LUT can share input) each lut can produce 2 outputs
    int adr_bits = ceil((double)log2(s));
    //int mux_tree_level = ceil((double)log2(mux_in)/2); //2 = log2(4)
    int mux_tree_level = ceil((double)adr_bits/2);//each level max 2 selects
    //int decoder_level = ceil((double)(adr_bits-1)/4); //fracturable
    int decoder_level = ceil((double)(adr_bits-1)/5);
    if (adr_bits == 1)
        decoder_level = 1;
    int mux_in = s;
    int dec_out = s;

    for (int i= 0; i<mux_tree_level ; i++)
    {
        if ((mux_in%4 == 1) ) //if only one input no lut needed
        {
            mux_in = mux_in/4+1;
            mux_lut += (mux_in-1); 
        }   
        else 
        {
            mux_in = ceil((double)mux_in/4);
            mux_lut += mux_in; 
        }
    }
    //cout << "decoder_level "<< decoder_level << endl;
    if (s<=2)
        dec_lut = 1;
    else 
    {
        for (int j =0; j<decoder_level; j++)
        {
            //dec_out = ceil((double)dec_out/2);//fracturable
            dec_out = ceil((double)dec_out);
            dec_lut += dec_out;
        }
    }
    //cout << "dec= " << dec_lut << " mux= " << mux_lut << endl;
    if (td)
        return 2*(dec_lut+mux_lut*lw);
    else 
        return dec_lut+mux_lut*lw;
}

///////////////////////////////////////////////////////////////// these to should be the same thing
int PhRAM::cost_to_map(RAM* r, int &parallel, int &serial, int &width) const 
{ //dt if RAM is dual port
    int d = r->depth;
    int w = r->width;

   return this->cost_to_map(d , w, (r->type == TrueDualPort) , parallel, serial, width);
}


int PhRAM::cost_to_map(int d , int w, bool td , int &parallel, int &serial, int &width) const 
{ //dt if RAM is dual port
    int max_width_p = max_width;
    if (td)
        max_width_p = max_width/2;
    int min_depth = max_size/max_width_p;

    //--------------------------------------------------------------------//
    //optimize for both
    int opt_w = 1; //min width 
    int opt_parallel = w;
    int opt_serial = (int)ceil((double)d/max_size); // base case if using in 1xSize(max/depth)
    if (opt_serial > MAX_SERIAL)
        return -1; // illigal //showing it can not be impelemented with this type
  
    // serial count will increase gradually as the width increase
    int opt_count_t = opt_parallel*opt_serial*ratio + added_logic_count(opt_serial,w,td); //w for max parralell , d / max/maxwidth for max serial
    for (int i = 1; i <= max_width_p ; i = i*2)
    {
        int c_parallel = (int)ceil((double)w/i); //how many in parallel
        int c_serial   = (int)ceil(d/((double)max_size/i)); // how many in series
        //cout << c_serial << " ";
        if (c_serial> MAX_SERIAL)
            break;
        //cout << "c_w =" << w/i << "c_d = " << d/(max_size/i) << " with i=" <<i <<   endl;
        //cout << "fot width " << i << " parallel=" << c_w << " serials=" << c_d << endl;
        int count_t = c_parallel*c_serial*ratio + added_logic_count(c_serial,w,td); // ! plus cost for dec and mux
        if (opt_count_t >= count_t) //serial is less costly than parallel so not using >=
        {
            opt_count_t = count_t;
            opt_serial = c_serial;
            opt_parallel = c_parallel;
            opt_w = i;
        }
    }
  
    serial = opt_serial;
    //cout << "opt=" << opt_serial << endl; 
    parallel = opt_parallel;
    width = opt_w;
    //cout << "optimium width is " << opt_i << endl; 
    //cout << "serial= " << opt_serial << " parallel "<< opt_parallel << endl; 
    //-----------------------------------------------------------------------//
    return opt_count_t;
}

///////////////////////////////////////////////////////////////// these to should be the same thing, different set of inputs

int LutRAM::cost_to_map(RAM* r, int &parallel, int &serial, int &width) const
{
    int d = r->depth;
    int w = r->width;

    this->cost_to_map(d, w, (r->type == TrueDualPort), parallel, serial, width);
}

int LutRAM::cost_to_map(int d, int w, bool td, int &parallel, int &serial, int &width) const
{
    int min_depth = max_size/max_width;
    int min_width = max_width/2; //10bit
    int max_depth = max_size/min_width;

    if (td)
        return -1; // illigal

    //if (support_LUTRAM == false )
        //return -1;

    int wider_parallel  = ceil((double)w/(max_width));
    int wider_serial    = ceil((double)d/(min_depth));

    int deeper_parallel = ceil((double)w/(min_width));
    int deeper_serial   = ceil((double)d/(max_depth));

    int opt_cost_t_wider  = wider_serial*wider_parallel*ratio   + added_logic_count(wider_serial,w,td); //32x20bit
    int opt_cost_t_deeper = deeper_serial*deeper_parallel*ratio + added_logic_count(deeper_serial,w,td); //64x10bit
    //int opt_count_t = min(opt_cost_t_wider, opt_cost_t_deeper);
    if ((opt_cost_t_wider <= opt_cost_t_deeper) && (wider_serial<MAX_SERIAL))
    {
        parallel = wider_parallel;
        serial   = wider_serial;
        width    = max_width;
        return opt_cost_t_wider;
    }
    else if  ((opt_cost_t_wider > opt_cost_t_deeper) && (deeper_serial<MAX_SERIAL))
    {
        parallel = deeper_parallel;
        serial   = deeper_serial;
        width    = min_width;
        return opt_cost_t_deeper;
    }
    else 
        return -1; //can not be impelemnted with LUTRAM

    /////////////////////////for debug
    // if (opt_cost_t_wider < opt_cost_t_deeper)
    //     cout << "serial= " << ceil((double)d/(min_depth)) << " parallel "<< ceil((double)w/(max_width)) << endl; 
    // else 
    //     cout << "serial= " << ceil((double)d/(max_depth)) << " parallel "<< ceil((double)w/(min_width)) << endl; 
    //////////////////////////
}

//===========================================================================//
Physical_RAM_Bank* Physical_RAM_Bank::copy()
{
    Physical_RAM_Bank* b = new Physical_RAM_Bank(id ,width, serial, parallel, mapped_logical_index, top_circuit,  type, logic_added);
    return b;
}

void Physical_RAM_Bank::modify(int s, int p, int w, int t , int lw, bool td)
{ 
    type = t;
    width = w;
    serial = s;
    parallel = p;
    //bool td  = (top_circuit->logical_RAMs[mapped_logical_index]->type == TrueDualPort)
    logic_added = added_logic_count(serial, lw, td); 
}
void Physical_RAM_Bank::modify(Physical_RAM_Bank* b)
{
    type = b->type;
    width = b->width;
    serial = b->serial;
    parallel = b->parallel; 
    logic_added = b->logic_added;
}

ostream& operator<<(ostream& out, const Physical_RAM_Bank* p){
    int index = p->mapped_logical_index;
    RAM* l = p->top_circuit->logical_RAMs[index];
    out << p->top_circuit->get_id() << " " << l->id 
        << " " << p->logic_added 
        << " LW " << l->width << " LD " << l->depth
        << " ID " << p->id << " S " << p->serial << " P " << p->parallel 
        << " Type " << (p->type+1)
        << " Mode " <<  RAMtype_name(l->type) //cause logical and physical are the same
        << " W " << p->width << " D " << phRAMs[p->type]->max_size/p->width << endl;
}

//===========================================================================//
int RAM::RAM_size() {
    //if depth > 256
    return width*depth;
}

bool RAM::equal(RAM* rc){
    if ((rc->width == width) && (rc->depth == depth))
        return true;
    return false; 
}  
//b is the target phRAM*2 and +1 if it is MOVEBIGGER
bool RAM::less(RAM* rc, int b){ //eqaul to higher position in sorted vector
    //cout << ">>>>>>>> comparing RAM "<< depth << " " << RAMtype_name(type) << " "<< width << " with this "<< rc   << endl;
    int c1[2],c2[2];
    int s1 , p1, w1;
    int s2 , p2, w2;
    int base = b/2;
    int traget_base;
    if (b%2 == 0 ) //MOVESMALLER
        traget_base = base-1;
    else //MOVEBIGGER 
        traget_base = base+1;

    c1[0] = phRAMs[base]->cost_to_map(    depth,    width,    (type==TrueDualPort), p1, s1, w1);
    c1[1] = phRAMs[traget_base]->cost_to_map(    depth,    width,    (type==TrueDualPort), p1, s1, w1);
    
    c2[0] = phRAMs[base]->cost_to_map(rc->depth,rc->width,(rc->type==TrueDualPort), p2, s2, w2);
    c2[1] = phRAMs[traget_base]->cost_to_map(rc->depth,rc->width,(rc->type==TrueDualPort), p2, s2, w2);

    //this should not happen
    /*if ((c1 != -1) && (c2 == -1) ) //pivot is illigal but ram is not => false //already there
        return true;
    if (c1 == -1)
        return false;*/
    //cout << c1[1] << "," << c1[0] << " and " << c2[1] << "," << c2[0] << endl;
    //cout << "is it less than = " << (c1[1]/c1[0] < c2[1]/c2[0]) << endl; 

    if (c1[0]==-1 || c1[1]==-1 || c2[0]==-1 || c2[1]==-1) 
    {
        cout << "error iligal values for comparing in RAM::less" << endl;
        throw Bad_Data();
    }

    return (c1[1]/c1[0] > c2[1]/c2[0]); //best case is at last element in the vector
    
    /*int c1[2];
    int c2[2];
    if (b==0) //LUTRAM
    {
        c1[1] = cost_map_LUTRAM(    depth,    width,    (type==TrueDualPort));//*phRAMs[LUTRAM].ratio;
        c2[1] = cost_map_LUTRAM(rc->depth,rc->width,(rc->type==TrueDualPort));//*phRAMs[LUTRAM].ratio;

        if ((c1[1] != -1) && (c2[1] == -1) ) //pivot is trudual but ram is not => false //already there
            return true;
        if (c1[1] == -1)
            return false;
        c2[0] = 1;
        c1[0] = 1;
    }
    else
    {
        c1[1] = cost_map_physical(    depth,    width,    (type==TrueDualPort),b); //*phRAMs[p].ratio;
        c2[1] = cost_map_physical(rc->depth,rc->width,(rc->type==TrueDualPort),b); //*phRAMs[p].ratio; 
        
        c1[0] = cost_map_physical(    depth,    width,    (type==TrueDualPort),b); 
        c2[0] = cost_map_physical(rc->depth,rc->width,(rc->type==TrueDualPort),b); 
    }
    //cout << "cost is " << c1 << "and " << c2 << endl;
    return (c1[1]/c1[0] < c2[1]/c2[0]);*/  
}

void RAM::map_lowest_individual_cost()
{
    int p , s, w;
    int c[MAX_PHRAM_TYPES_COUNT];
    int min_cost;
    int min_cost_PRAM;
    bool init = false;
    for (int r = 0; r < PHRAM_TYPES_COUNT ; r++)
    {
        c[r] = phRAMs[r]->cost_to_map(depth, width, (type == TrueDualPort), p , s, w);
        if (c[r] != -1 ) //if -1 should not use
        {
            if (init == false)
            {
                min_cost = c[r];
                min_cost_PRAM = r; 
                init = true;  
            }
            else if (c[r] <= min_cost) //= because it means less memory needed
            {
                min_cost = c[r];
                min_cost_PRAM = r;
            }
        }
    }
    //cout << this;
    //cout << c[0] << "  " << c[1]<< "  " << c[2] << endl;
    //cout << "min is >>" << c[min_cost_PRAM] << " index =" << min_cost_PRAM << endl;
    if (init == false)
    {
        cout << "error no choice for mapping, all having parallel greater than 16" << endl;
        throw Bad_Data();
    }
    //updating values of p,s,w
    phRAMs[min_cost_PRAM]->cost_to_map(depth, width, (type == TrueDualPort), p , s, w);
    Physical_RAM_Bank* new_physical_RAM = new Physical_RAM_Bank(top_circuit->physical_RAMs.size() ,w, s, p, id, top_circuit,  min_cost_PRAM, added_logic_count(s, width, (type == TrueDualPort)));
    
    top_circuit->physical_RAMs.push_back(new_physical_RAM);
    if (min_cost_PRAM > 0)
    {
        if (c[min_cost_PRAM-1] != -1)
        {
            top_circuit->sorted_logical_RAMs[min_cost_PRAM][MOVESMALLER].push_back(id);
            //cout << "adding this one to samller, id= " << id << endl;
        }
    }
    if (min_cost_PRAM < (PHRAM_TYPES_COUNT-1))
    {
        if (c[min_cost_PRAM+1] != -1)
        {
            top_circuit->sorted_logical_RAMs[min_cost_PRAM][MOVEBIGGER].push_back(id);
            //cout << "adding this one to bigger, id= " << id << endl;
        }
    }
    mapped_RAM = new_physical_RAM;
    mapped = true;
}

// string RAM::checker_output()
// {
//     stringstream ss;
//     ss << " LW " << width << " LD " << depth <<
//}

ostream& operator<<(ostream& out, const RAM& r) {
    out << "id: "<<r.id << "\t" << RAMtype_name(r.type) << "\t depth: " << r.depth << " width: "<< r.width << endl;
    return out;
}

ostream& operator<<(ostream& out, const RAM* r) {
    out << "id: "<<r->id << "\t" << RAMtype_name(r->type) << "\t depth: " << r->depth << " width: "<< r->width << endl;
    return out;
}
//=====================================================================//
void Circuit::add_logical_RAM(RAM* new_RAM) {
    logical_RAMs.push_back(new_RAM);
    for (int r = 0; r < PHRAM_TYPES_COUNT ; r++)
    {
        //cout << new_RAM->id << " into " << r << endl;
        //sorted_logical_RAMs[r].push_back(new_RAM->id); // if all want all the RAMs uncomment
    }
    if ( new_RAM->type == TrueDualPort)
        unmapped_tdual.push_back(new_RAM->id);
}

ostream& operator<<(ostream& out, const Circuit& c) {
    out << "Circuit # " << c.circuit_id << " : " << endl;
    out << '\t' << "logic block count = " << c.logic_block_count <<endl;
    for (int i = 0; i < c.logical_RAMs.size(); i++)
        out << '\t' << c.logical_RAMs[i];
    return out;
}
ostream& operator<<(ostream& out, const Circuit* c) {
    out << "Circuit # " << c->circuit_id << " : " << endl;
    out << '\t' << "logic block count = " << c->logic_block_count <<endl;
    for (int i = 0; i < c->logical_RAMs.size(); i++)
        out << '\t' << c->logical_RAMs[i];
    for (int j=0; j<PHRAM_TYPES_COUNT; j++)
    {
        out << "@ sorted for RAM TYPE " << PHRAMtype_name(phRAMs[j]->type) << endl;
        out << "\tBIGGER" << endl;
        for (int i = 0; i < c->sorted_logical_RAMs[j][MOVEBIGGER].size(); i++)
            out << '\t' << c->logical_RAMs[c->sorted_logical_RAMs[j][MOVEBIGGER][i]];
        out << "\tSMALLER" << endl;
        for (int i = 0; i < c->sorted_logical_RAMs[j][MOVESMALLER].size(); i++)
            out << '\t' << c->logical_RAMs[c->sorted_logical_RAMs[j][MOVESMALLER][i]];
    }
    return out;
}

void Circuit::update_available_resources(){
    for (int j =  0; j < PHRAM_TYPES_COUNT; j++)
    {
        if (LutRAM* lr = dynamic_cast<LutRAM*>(phRAMs[j])) // then it is LUTRAM
            available_resources[j] = 0; //adding LUTRAM without penalty
        else
            available_resources[j] = (int)(logic_block_count/phRAMs[j]->ratio);
    }
    //cout << available_resources[LUTRAM] << " M8K =" <<available_resources[M8K]  << " M128="<< available_resources[M128K] << endl;
}


void Circuit::update_required_resources()
{
    for (int j =0; j<= PHRAM_TYPES_COUNT ; j++)
        required_resources[j] = 0; //also update for logical 

    for (int p =0; p< physical_RAMs.size(); p++ )
    {
        Physical_RAM_Bank* ram_bank = physical_RAMs[p];
        int t = ram_bank->type;
        required_resources[t] += (ram_bank->serial)*(ram_bank->parallel);
        //required_area[t] = 
        required_resources[PHRAM_TYPES_COUNT] += ram_bank->logic_added; 
    }
    required_resources[PHRAM_TYPES_COUNT]  = ceill((double)required_resources[PHRAM_TYPES_COUNT]/N);
    if (ENABLE_COUT)
    {
        cout << "required: TYPE0= " <<required_resources[0] << " TYPE1 =" <<required_resources[1]  << " TYPE2="<< required_resources[2] << " and added logic= " << required_resources[PHRAM_TYPES_COUNT]  << endl;
        cout << "available: TYPE0 = "<< available_resources[0] << " TYPE1 =" <<available_resources[1]  << " TYPE2="<< available_resources[2] << " and logic = " <<  logic_block_count <<  endl;
    }
}

int Circuit::update_max_resources(int& index){
    this->update_required_resources();
    int max_area = 0;
    int max_area_index = 0;
    int a[MAX_PHRAM_TYPES_COUNT];
    int logical_LUT = logic_block_count + required_resources[PHRAM_TYPES_COUNT];
    for ( int k =0 ; k<PHRAM_TYPES_COUNT ; k++)
    {
        a[k] = required_resources[k]*phRAMs[k]->ratio;
        if (LutRAM* lr = dynamic_cast<LutRAM*>(phRAMs[k]))
        {
            if ( a[k] < logical_LUT+ (a[k]/phRAMs[k]->ratio) )  //logic_LUT+ (a[LUTRAM]/phRAMs[LUTRAM].ratio) ==> total lut needed
                a[k] =  logical_LUT+ (a[k]/phRAMs[k]->ratio);
            else //more LUTRAM 
                a[k] =  a[k] - logical_LUT;// (a[LUTRAM]-logic_block_count) = empty blocks
        }

        if (a[k] > max_area)
        {
            max_area_index = k;
            max_area = a[k];
        }
    }
    
    for (int i = 0 ; i<PHRAM_TYPES_COUNT; i++ )
        max_available_resources[i] = max_area/phRAMs[i]->ratio;

    geo_area = 0;
    for (int i = 0 ; i<PHRAM_TYPES_COUNT; i++ )
        geo_area += max_available_resources[i]*phRAMs[i]->area;


    
    index = max_area_index;
    return max_area;
    //return geo_area;
}

void Circuit::update_area()
{
    area = this->update_max_resources(limiting_rsc);
    if (ENABLE_COUT)
        cout << "area = " << area << " (of numer of LUTs) and limiting is "<< PHRAMtype_name(limiting_rsc) << endl;

}

bool my_insert(vector<int> &slr, int index, int t, int big_small, Circuit* c, RAM* r)
{
    int p,s,w,j;
    int target = 2*t;
    int target2 = t-1;
    if (big_small == MOVEBIGGER)
    {
        target++;
        target2 = t+1;
    }

    if (t >= PHRAM_TYPES_COUNT || t<0 || 
        target<0 || target>= PHRAM_TYPES_COUNT||
        target2<0 || target2>= PHRAM_TYPES_COUNT )
        return false;

    int cost = phRAMs[target2]->cost_to_map(r, p, s, w);
    if (cost == -1)
        return false; 

    for ( j =0; j < slr.size(); j++)
    {
        if (c->logical_RAMs[slr[j]]->less(r, target)) //!!!!
            break;
    }
    //  c->sorted_logical_RAMs[t][big_small].push_back(0);
    //  for (int i = c->sorted_logical_RAMs[t][big_small].size() ; i > j ; i--)
    //  {
    //     c->sorted_logical_RAMs[t][big_small][i] = c->sorted_logical_RAMs[t][big_small][i-1];
    //  }
    // c->sorted_logical_RAMs[t][big_small][j] = index;
    c->sorted_logical_RAMs[t][big_small].insert(c->sorted_logical_RAMs[t][big_small].begin()+j, index); //insert vlaue i
}
bool my_erase(vector<int> &slr, int i)
{
    for ( int j =0; j < slr.size(); j++)
    {
        if (slr[j] == i)
        {
            slr.erase(slr.begin()+j);
        }

    } 
}

bool Circuit::modify_mapping(int t, int Small_BIG) //t is source target is t+1
{
    int p,s,w;
    int target;
    int Not_;

    if ((t < 0) || (t >= PHRAM_TYPES_COUNT))
    {
        cout << "Error in map modify t is" << t << endl;
        throw Bad_Data();
    }

    if ( sorted_logical_RAMs[t][Small_BIG].empty())
        return false;

    if (Small_BIG == MOVEBIGGER)
    {
        target = t+1;
        Not_ = MOVESMALLER;
    }
    else if (Small_BIG == MOVESMALLER)
    {
        target = t-1;
        Not_ = MOVEBIGGER;
    }
    else 
        throw Bad_Data();


    //if it is then no benefit to do it
    //this should never happen because then it is not limiting resource
    int index = sorted_logical_RAMs[t][Small_BIG].back(); 
    //phRAMs[t].cost_to_map(logical_rams[index], p_prev, s_prev, w_prev);
    phRAMs[target]->cost_to_map(logical_RAMs[index], p, s, w);
    if (ENABLE_COUT)
        cout << "RAM chosen to move from " << PHRAMtype_name(t) <<  " to "<< PHRAMtype_name(target) << ": " << logical_RAMs[index] ;
    Physical_RAM_Bank* b = logical_RAMs[index]->mapped_RAM->copy();
    //cout << "copied" << b;
    //cout << ">>>>>>logical width will be "<< w << "and parallel =" << p << endl;
    logical_RAMs[index]->mapped_RAM->modify(s,p,w, target,logical_RAMs[index]->width , (logical_RAMs[index]->type == TrueDualPort));
    int prev_area = area;
    update_area();
    if (area < prev_area) //good, keep it
    {

        sorted_logical_RAMs[t][Small_BIG].pop_back();
        my_erase(sorted_logical_RAMs[t][Not_], index);
        sorted_logical_RAMs[target][Not_].push_back(index);
        //my_insert(sorted_logical_RAMs[target][Not_], index, target, this, logical_RAMs[index]);
        my_insert(sorted_logical_RAMs[target][Small_BIG], index, target, Small_BIG, this, logical_RAMs[index]);
        delete b;
        if (ENABLE_COUT)
            cout << "improv made" << endl;
        return true;
    }
    else // reverse what has done
    {
        logical_RAMs[index]->mapped_RAM->modify(b);
        delete b;
        update_area();
        if (area != prev_area)
        {
            cout << "Error in function move_bigger" << endl;
            throw Bad_Data();
        }
        return false;
    }
}

void Circuit::optmize_area()
{
    this->update_area(); //it will update max_resoucers_also
    int move_big_small, move_big_small_not;
    while (1)
    {
        int index = limiting_rsc;
        if (required_resources[index] <= available_resources[index])
            break;
        
        if (index != PHRAM_TYPES_COUNT-1) //first choose to try for bigger available one if not go for samller 
        {
            if (required_resources[index+1] < available_resources[index+1])
            {
                move_big_small = MOVEBIGGER;
                move_big_small_not = MOVESMALLER;
            }
            else 
            {
                move_big_small = MOVESMALLER;
                move_big_small_not = MOVEBIGGER;
            }
        }
        else 
        {
            move_big_small = MOVESMALLER;
            move_big_small_not = MOVEBIGGER;
        }

        if (modify_mapping(index,move_big_small))
            continue;
        else if (modify_mapping(index, move_big_small_not))
            continue;
        else //no furture improv possible 
            break;
    } 
}

Circuit::~Circuit(){
    for (int i = 0; i < logical_RAMs.size(); i++)
        delete logical_RAMs[i];
    for (int i = 0; i < physical_RAMs.size(); i++)
        delete physical_RAMs[i];
}
//================================================================================//

vector<Circuit*> read_LBC_file(string file_name){
    vector<Circuit*> circuits;
    int id, LBcount;
    string line;
    ifstream LBC_file (file_name.c_str());
    if (LBC_file.is_open())
    {
        getline (LBC_file,line); // first line for explanantion
        while ( getline (LBC_file,line) )
        {
            istringstream l(line);
            l>>id;
            l>>LBcount;
            Circuit* new_circuit = new Circuit(id,LBcount);
            circuits.push_back(new_circuit);
            circuits.back()->update_available_resources();

        }
        LBC_file.close();
        return circuits;
    }
    else 
    {
        cout << "unable to open file " << file_name << endl ;
        throw Bad_Data();
    }
}

bool read_LR_file(string file_name, vector<Circuit*>& circuits){
    int c_id, r_id, W, D;
    string t; 
    string line;
    ifstream LR_file (file_name.c_str());
    if (LR_file.is_open())
    {
        getline (LR_file,line); // first line for number of circuit counts
        getline (LR_file,line); // first line for explanantion
        while ( getline (LR_file,line) )
        {
            istringstream l(line);
            l>>c_id;
            l>>r_id;
            l>>t;
            l>>D;
            l>>W;
            RAM* new_RAM = new RAM(r_id, t, D, W, circuits[c_id]);
            //cout << "circuit" << c_id << "-> ";
            //cout << new_RAM ;
            circuits[c_id]->add_logical_RAM(new_RAM);
        }
        LR_file.close();
    }
    else 
    {
        cout << "unable to open file " << file_name << endl;
        throw Bad_Data();
    }
}

bool write_PR_file(string file_name, vector<Circuit*> circuits){
    ofstream PR_file (file_name.c_str(), ios::trunc);
    if (PR_file.is_open())
    {
        for ( int j =0; j < circuits.size(); j++)
        {
            for (int p =0; p< circuits[j]->physical_RAMs.size(); p++ )
                PR_file << circuits[j]->physical_RAMs[p];
        }
        PR_file.close();
    }
    else 
    {
        cout << "unable to open file " << file_name << endl;
        throw Bad_Data();
    }
}

// A utility function to swap two RAMs
void swap_RAM(vector<RAM*> &RAMs, int i, int j)
{
    RAM* tmp = RAMs[j];
    RAMs[j] = RAMs[i];
    RAMs[i] = tmp;
}
// A utility function to swap two RAMs
void swap(vector<int> &RAMs, int i, int j)
{
    int tmp = RAMs[j];
    RAMs[j] = RAMs[i];
    RAMs[i] = tmp;
}

int partition (vector<RAM*> LRAMs,vector<int> &sorted_RAMs, int left, int right, int base, int &repetition)
{
    /* This function takes random RAM as pivot, places the pivot RAM at its correct position in sorted
    vector, and places all smaller (smaller than pivot in total size in bits) to left of pivot and
    all greater in total size in bits RAMs to right of pivot */
    //randomly chosing pivot because files are kind of sorted and for optimization, balanced devision is needed
    //for faster sort because of often repetition sorting out them 
    
    int rand_index = left + (rand() % (right - left+1));
    swap(sorted_RAMs , right ,rand_index);
    int pivot = sorted_RAMs[right];    // pivot
    RAM* pivot_RAM = LRAMs[pivot];
    int i = (left - 1);  // where exactly pivot happens in data// Index of smaller element
    int k = (left - 1);
 
    for (int j = left; j < right; j++)
    {
        // If current element is smaller than or
        // equal to pivot
        if ( LRAMs[sorted_RAMs[j]]->equal(pivot_RAM) ) //repetition //if it costs less
        {
            i++;    // increment index of smaller element
            swap(sorted_RAMs, i, j);
            k++;
            swap(sorted_RAMs, i, k);
        }
        else if (LRAMs[sorted_RAMs[j]]->less(pivot_RAM, base)) //(RAMs[j]->RAM_size() > pivot) //shoud sort downer
        {
            //cout << "&&&&it is less than " << endl;
            i++;    // increment index of smaller element
            swap(sorted_RAMs, i, j);
        }
    }
    swap(sorted_RAMs , (i + 1) , right); //puting pivot in its place 
    for (int e = left  ; e <= k ; e++)
    {
        swap(sorted_RAMs,e, (i - (e - left)));
    }
    //cout << "pivot is : " << pivot << " i is: " << i << " rand: " << rand_index << " range " << left  <<"-" << right << " RAM: " << LRAMs[i + 1] ;
    repetition = k - (left-1);
    return (i + 1);
}
 
/* The main function that implements QuickSort*/
void quickSort(vector<RAM*> LRAMs, vector<int> &sorted_RAMs, int left, int right, int base)
{
    int rp;
    if (left < right)
    {
        /* pi is partitioning index */
        int pi = partition(LRAMs,sorted_RAMs, left, right, base, rp);
 
        // Separately sort elements before
        // partition and after partition
        quickSort(LRAMs,sorted_RAMs, left  , pi - 1 - rp,base);
        quickSort(LRAMs,sorted_RAMs, pi + 1, right      ,base);
    }
}

bool sort_circuit_RAM (Circuit* c)
{
    for (int r = 0; r < PHRAM_TYPES_COUNT ; r++)
        {       
            if (r>0) // nothing less than LUTRAM
                    quickSort(c->logical_RAMs,c->sorted_logical_RAMs[r][MOVESMALLER],0,c->sorted_logical_RAMs[r][MOVESMALLER].size()-1, 2*r);
            if (r<(PHRAM_TYPES_COUNT-1))
                    quickSort(c->logical_RAMs,c->sorted_logical_RAMs[r][MOVEBIGGER],0,c->sorted_logical_RAMs[r][MOVEBIGGER].size()-1, 2*r+1);
        }

}

bool read_config(int argc, char* argv[])
{
    int it = 1; //skip the exe name 
    int r;
    LutRAM* new_LUT_RAM; 
    PhRAM* new_physical_RAM;
    while(it < argc - 3)
    {
        if (argv[it][0] == '-')
        {
            switch (argv[it][1]){
                case'd':
                        for (int p=0; p<MAX_PHRAM_TYPES_COUNT; p++)
                        {
                            STRATIX = true;
                            if (p==0)
                            {
                                new_physical_RAM = new LutRAM();
                                new_physical_RAM->area = (LB_AREA_LUT+LB_AREA)/LB_COUNT_LUTRAM;
                            }
                            else 
                                new_physical_RAM = new PhRAM(p);
                            phRAMs[p]= new_physical_RAM;
                            //phRAMs.push_back(new_physical_RAM);
                        }
                        it++;
                        PHRAM_TYPES_COUNT = MAX_PHRAM_TYPES_COUNT;
                        break;
                case 'l':
                        it = it+2; //1:1
                        r = atoi(argv[it]); 
                        it++;
                        new_LUT_RAM = new LutRAM();
                        new_LUT_RAM->area = (LB_AREA_LUT+LB_AREA)/r;
                        // if (argv[it][0] == 'L')
                        //     new_LUT_RAM->support_LUTRAM = false;
                        // else if (argv[it][0] == 'R')
                        //     new_LUT_RAM->support_LUTRAM = true;
                        // else
                        //     throw Bad_Data();
                        new_LUT_RAM->ratio = r;
                        //phRAMs.push_back(new_LUT_RAM);
                        phRAMs[PHRAM_TYPES_COUNT] = new_LUT_RAM;
                        PHRAM_TYPES_COUNT++;
                        break;
                case 'b':
                        it++;
                        new_physical_RAM = new PhRAM(atoi(argv[it]),atoi(argv[it+1]),atoi(argv[it+2]),PHRAM_TYPES_COUNT);
                        //phRAMs.push_back(new_physical_RAM);
                        phRAMs[PHRAM_TYPES_COUNT] = new_physical_RAM;
                        it = it+4;
                        PHRAM_TYPES_COUNT++;

                        break;
                default :
                        cout << "this option not supported, try d,b,l" << endl;
                        throw Bad_Data();
            }
        }
        else 
        {
            cout << "check commandline inputs" << endl;
            throw Bad_Data();
        }
    }
    //RAM_TYPES_COUNT = phRAMs.size();

}

int main(int argc, char* argv[]) {
    clock_t start_t, end_t;
    start_t = clock();
    string file_names[3];
    if (argc<5 )
    {
        cout << "check commandline!\n";
        throw Bad_Data();
    }
    file_names[2] = argv[argc-1];
    file_names[1] = argv[argc-2];
    file_names[0] = argv[argc-3];
    //cout << "file_names are " << file_names[0] << file_names[1] << file_names[2]<< endl;
    srand(time(0)); //for quick sort
    //LutRAM(),PhRAM(M8K),PhRAM(M128K)
    //-l 1 1 R -b 8192 32 10 1 -b 131072 128 300 1 logical_rams.txt logic_block_count.txt physical_rams.txt 

    read_config(argc,argv);

    //cout<< phRAMs[0]->max_size << " width= " << phRAMs[0]->max_width << " and type = " << phRAMs[0]->type << " and ratio = " << phRAMs[0]->ratio  << " and area " << phRAMs[0]->area << endl;
    //cout<< phRAMs[1]->max_size << " width= " << phRAMs[1]->max_width << " and type = " << phRAMs[1]->type << " and ratio = " << phRAMs[1]->ratio  << " and area " << phRAMs[1]->area << endl;
    //cout<< phRAMs[2]->max_size << " width= " << phRAMs[2]->max_width << " and type = " << phRAMs[2]->type << " and ratio = " << phRAMs[2]->ratio  << " and area " << phRAMs[2]->area << endl;
    //cout << cost_map_physical(40900, 14, 1, M8K);
    vector<Circuit*> circuits = read_LBC_file(file_names[1]);
    read_LR_file(file_names[0],circuits);
    //for (int j =0; j<circuit_count; j++)
    //cost_map_LUTRAM(9,4,0);

    //phRAMs[LUTRAM].cost_to_map(1024,16,0, p, s);
    //cout << "LUTRAM " << s << " and " << p << endl;
    //phRAMs[M8K].cost_to_map(1024,16,0, p, s);
    //cout << "M8K " << s << " and " << p << endl;
    //phRAMs[M128K].cost_to_map(1024,16,0, p, s);
    //cout << "M128K " << s << " and " << p << endl;
    
    //for (int j =0; j<circuits.size(); j++)
    //        cout << circuits[j] ;
    
    double geometric_area = 1;
    //int j=63;
    for (int j =0; j<circuits.size(); j++)
    {
        for (int k = 0; k < circuits[j]->logical_RAMs.size(); k++ )
        {
            circuits[j]->logical_RAMs[k]->map_lowest_individual_cost();
        }
        circuits[j]->update_area();
        sort_circuit_RAM(circuits[j]);
        if (ENABLE_COUT)
            cout << circuits[j];
        circuits[j]->optmize_area();
        write_PR_file(file_names[2],circuits);
        geometric_area = geometric_area*(circuits[j]->get_geo_area()/(pow(10,8)));

        //cout << circuits[j];
        //for (int p =0; p< circuits[j]->physical_RAMs.size(); p++ )
        //    cout << circuits[j]->physical_RAMs[p];
        //cout << circuits[j];
    }

    cout << "area is " << pow(geometric_area, 1.0/circuits.size())*pow(10,8) << endl;


    //cout << added_logic_count(100, 1) << endl;


    //============================================================
    // for (int j =62; j<64; j++)
    // {
    //     for (int r = 0; r < PHRAM_TYPES_COUNT ; r++)
    //     {
    //         quickSort(circuits[j]->logical_RAMs,circuits[j]->sorted_logical_RAMs[r],0,circuits[j]->logical_RAMs.size()-1, r);
    //         for (int k = 0; k < circuits[j]->sorted_logical_RAMs[r].size(); k++ )
    //             cout << circuits[j]->sorted_logical_RAMs[r][k] << " ";
    //         cout << endl;
    //     }
    //     cout << circuits[j] << "**************************************\n";
    // }
    end_t = clock();
    double total_t = (double)(end_t - start_t);

    cout << "Total clks taken by CPU:" << total_t << endl;
}