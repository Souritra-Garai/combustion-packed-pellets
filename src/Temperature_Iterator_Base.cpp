#include "Temperature_Iterator_Base.hpp"

// Function definitions for class Temperature Iterator Base

Temperature_Iterator_Base::Temperature_Iterator_Base(unsigned int n)
{
    N = n;
    T_ign = 0;

    reset_iterators();
}

void Temperature_Iterator_Base::reset_iterators()
{ 
    i = 0;
    reset_other_iterators();
}

bool Temperature_Iterator_Base::iterators_in_range() { return i < N; }

bool Temperature_Iterator_Base::in_reaction_zone(float T)
{
    return (T > T_ign && iterators_in_range());
}

bool Temperature_Iterator_Base::in_post_combustion_zone(float eta)
{
    return (eta > 0.99999 && iterators_in_range());
}

void Temperature_Iterator_Base::increment_iterators()
{ 
    i++;
    increment_other_iterators();
}

void Temperature_Iterator_Base::increment_iterators(
    vector<float>::iterator i1  )
{
    increment_iterators();
    i1++;
}

void Temperature_Iterator_Base::increment_iterators(
    vector<float>::iterator i1,
    vector<float>::iterator i2  )
{
    increment_iterators();
    i1++;
    i2++;
}

void Temperature_Iterator_Base::increment_iterators(
    vector<float>::iterator i1,
    vector<float>::iterator i2,
    vector<float>::iterator i3  )
{
    increment_iterators();
    i1++;
    i2++;
    i3++;
}

void Temperature_Iterator_Base::increment_iterators(
    vector<float>::iterator i1,
    vector<float>::iterator i2,
    vector<float>::iterator i3,
    vector<float>::iterator i4  )
{
    increment_iterators();
    i1++;
    i2++;
    i3++;
    i4++;
}

void Temperature_Iterator_Base::reset_other_iterators() {;}

void Temperature_Iterator_Base::increment_other_iterators() {;}

void Temperature_Iterator_Base::set_ignition_temperature(float T)
{
    T_ign = T;
}
