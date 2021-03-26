#ifndef __TEMP_ITER_BASE__
#define __TEMP_ITER_BASE__

#include <vector>

using namespace std;

class Temperature_Iterator_Base
{
    protected:

        unsigned int N;
        
        float T_ign;
        unsigned int i;

        void reset_iterators();
        
        bool iterators_in_range();
        bool in_reaction_zone(float T);
        bool in_post_combustion_zone(float eta);

        void increment_iterators();
        void increment_iterators(vector<float>::iterator);
        void increment_iterators(vector<float>::iterator, vector<float>::iterator);
        void increment_iterators(vector<float>::iterator, vector<float>::iterator, vector<float>::iterator);
        void increment_iterators(vector<float>::iterator, vector<float>::iterator, vector<float>::iterator, vector<float>::iterator);

        virtual void reset_other_iterators();
        virtual void increment_other_iterators();

    public:

        Temperature_Iterator_Base(unsigned int n);
        void set_ignition_temperature(float T_Ignition);
};

#endif