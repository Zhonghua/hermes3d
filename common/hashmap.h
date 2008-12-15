#ifndef _MAP_H_
#define _MAP_H_

#include <assert.h>

/// \file map.h
/// \brief Template for a hash map.

#include <Judy.h>

#ifndef JUDY_INVALID
#define JUDY_INVALID					((Word_t) -1)
#endif

/// \class Map
/// \brief Implementation of dynamic array, using an array-of-bytes of Length as an Index and a word as a Value.
/// C++ encapsulation of JudyHS functions.
/// Provides functionality of a hash map.
/// Items are classes of TYPE
template<class TYPE>
class Map {
protected:
	void *judy_hs;
	void *judy_l;

public:
	Map();
	virtual ~Map();

	int count();

	bool is_empty();

	bool lookup(char *key, int length, TYPE &item);

	void set(char *key, int length, TYPE item);

	bool remove(char *key, int length);

	void remove_all();

	// Iterators

	/// Get the first index that is present and is equal to or greater than the passed \c idx.
	/// Typically used to begin an iteration over all indices present in the array.
	/// \param[in] idx Optional, default value \c 0 (finds the first present index).
	/// \return
	/// 	\li First index present in the array that is equal or greater than the passed \c idx (if found), 
	/// 	\li \c JUDY_INVALID (if not found).
	Word_t first(Word_t idx = 0) const;

	/// Get the first index that is present and is greater than the passed \c idx.
	/// Typically used to continue an iteration over all indices present in the array.
	/// \param[in] idx Index whose succesor we want to find. Optional, default value \c 0.
	/// \return
	/// 	\li First idx present in the array that is greater than the passed \c idx (if found), 
	/// 	\li \c JUDY_INVALID (if not found).
	Word_t next(Word_t idx = 0) const;
	
	/// Get the last index present in the array that is equal to or less than the passed \c idx.
	/// Typically used to begin a reverse iteration over all indices present in the array.
	/// \param[in] idx Optional, default value <c>(Word_t) -1</c> (finds the last index present in the array).
	/// \return
	///		\li Last index present in the array that is equal or less than the passed \c idx (if found), 
	/// 	\li \c JUDY_INVALID (if not found).
	Word_t last(Word_t idx = (Word_t) -1) const;
	
	/// Get the last index present in the array that is less than the passed \c idx.
	/// Typically used to continue a reverse iteration over all indices present in the array.
	/// \param[in] idx Index whose predecessor we want to find. Optional, default value <c>(Word_t) -1</c>.
	/// \return
	/// 	\li Last index present in the array that is less than the passed \c idx (if found), 
	/// 	\li \c JUDY_INVALID (if not found).
	Word_t prev(Word_t idx = (Word_t) -1) const;
};

// implementation

template<class TYPE>
Map<TYPE>::Map() {
	judy = NULL;
};


#endif
