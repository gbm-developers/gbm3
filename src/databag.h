//------------------------------------------------------------------------------
//
//  File:       databag.h
//
//  Description:   header file for bagged data
//
//  Owner:      James Hickey
//
//  History:    07/06/2016  jhickey created.
//
//------------------------------------------------------------------------------

#ifndef DATABAG_H
#define DATABAG_H

//------------------------------
// Includes
//------------------------------
#include "datadistparams.h"

//------------------------------
// class definition
//------------------------------
class Bag {
public:
	//----------------------
	// Public Constructors
	//----------------------
	Bag(const DataDistParams& dataparams):
		is_in_bag_(dataparams.num_trainrows, false) {};

	//----------------------
	// Public Destructors
	//----------------------
	~Bag() {};
	//----------------------
	// Public Functions
	//----------------------
	bool get_element(long index) const { return is_in_bag_[index]; }
	void set_element(long index) { is_in_bag_[index] = 1; };
	void clear() { is_in_bag_.assign(is_in_bag_.size(), 0); };

private:
	//----------------------
	// Private Variables
	//----------------------
	std::vector<int> is_in_bag_;

};
#endif // DATABAG_H
