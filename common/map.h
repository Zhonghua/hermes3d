#ifndef _MAP_H_
#define _MAP_H_

#include <assert.h>

#include <Judy.h>

#ifndef INVALID_IDX
#define INVALID_IDX					((Word_t) -1)
#endif

/// \class Map
/// \brief Implementation of dynamic array, using an array-of-bytes of Length as an Index and a word as a Value.
/// Provides functionality of a hash map.
/// Items are classes of TYPE
template<class KEY, class TYPE>
class Map {
protected:
	void *judy_hs;
	void *judy_l;

public:
	Map();
	virtual ~Map();

	// Get the number of values in the Map
	// @return number of (key, item) pairs
	Word_t count() const;

	/// Test if the Map is empty.
	/// @return true if the Map has no items, otherwise false
	bool is_empty() const;

	bool lookup(const KEY &key, TYPE &item) const;

	Word_t get_idx(const KEY &key) const;

	/// Get a value as \c iter position
	/// \param[in] iter Iterator obtained by first(), next(), last() or prev().
	/// \return item at position \c iter
	TYPE get(Word_t iter) const;

	TYPE operator[](int idx) const;

	/// Add a new (key, item) pair into the Map
	/// \param[in] key Pointer to the array-of-bytes.
	/// \param[in] length Length of \c key
	/// \param[in] item Item to insert
	bool set(const KEY &key, TYPE item);

	/// Delete an item with key \c key from the Map.
	/// \param[in] key Pointer to the array-of-bytes.
	/// \param[in] length Length of \c key
	bool remove(const KEY &key);

	/// Remove all items from the array
	/// Does NOT destruct items in the array
	void remove_all();

	// Iterators

	// NOTE: during the iteration, keys of items are not available

	/// Get the first index that is present and is equal to or greater than the passed \c idx.
	/// Typically used to begin an iteration over all indices present in the array.
	/// \param[in] idx Optional, default value \c 0 (finds the first present index).
	/// \return
	/// 	\li First index present in the array that is equal or greater than the passed \c idx (if found),
	/// 	\li \c INVALID_IDX (if not found).
	Word_t first() const;

	/// Get the first index that is present and is greater than the passed \c idx.
	/// Typically used to continue an iteration over all indices present in the array.
	/// \param[in] idx Index whose succesor we want to find. Optional, default value \c 0.
	/// \return
	/// 	\li First idx present in the array that is greater than the passed \c idx (if found),
	/// 	\li \c INVALID_IDX (if not found).
	Word_t next(Word_t idx) const;

	/// Get the last index present in the array that is equal to or less than the passed \c idx.
	/// Typically used to begin a reverse iteration over all indices present in the array.
	/// \param[in] idx Optional, default value <c>(Word_t) -1</c> (finds the last index present in the array).
	/// \return
	///		\li Last index present in the array that is equal or less than the passed \c idx (if found),
	/// 	\li \c INVALID_IDX (if not found).
	Word_t last() const;

	/// Get the last index present in the array that is less than the passed \c idx.
	/// Typically used to continue a reverse iteration over all indices present in the array.
	/// \param[in] idx Index whose predecessor we want to find. Optional, default value <c>(Word_t) -1</c>.
	/// \return
	/// 	\li Last index present in the array that is less than the passed \c idx (if found),
	/// 	\li \c INVALID_IDX (if not found).
	Word_t prev(Word_t idx) const;

protected:
	void free_item(Word_t idx);
};

// implementation

template<class KEY, class TYPE>
Map<KEY, TYPE>::Map() {
	judy_hs = NULL;
	judy_l = NULL;
};

template<class KEY, class TYPE>
Map<KEY, TYPE>::~Map() {
	remove_all();
};

template<class KEY, class TYPE>
Word_t Map<KEY, TYPE>::count() const {
	Word_t count;
	JLC(count, judy_l, 0, -1);
	return count;
}

template<class KEY, class TYPE>
bool Map<KEY, TYPE>::is_empty() const {
	return count() == 0;
}

template<class KEY, class TYPE>
bool Map<KEY, TYPE>::lookup(const KEY &key, TYPE &item) const {
	void *pval = NULL;
	// check if the key exists
	JHSG(pval, judy_hs, (char *) &key, sizeof(KEY));
	if (pval == NULL)
		return false;
	else {
		// get associated value
		Word_t idx = *(Word_t *) pval;
		JLG(pval, judy_l, idx);
		if (pval == NULL)
			return false;
		else {
			item = **((TYPE **) pval);
			return true;
		}
	}
}

template<class KEY, class TYPE>
Word_t Map<KEY, TYPE>::get_idx(const KEY &key) const {
	void *pval = NULL;
	// check if the key exists
	JHSG(pval, judy_hs, (char *) &key, sizeof(KEY));
	if (pval == NULL)
		return INVALID_IDX;
	else
		// get associated value
		return *(Word_t *) pval;
}

template<class KEY, class TYPE>
TYPE Map<KEY, TYPE>::get(Word_t iter) const {
	void *pval = NULL;
	JLG(pval, judy_l, iter);
	assert(pval != NULL);
	return **((TYPE **) pval);
}

template<class KEY, class TYPE>
TYPE Map<KEY, TYPE>::operator[](int idx) const {
	return get(idx);
}

template<class KEY, class TYPE>
bool Map<KEY, TYPE>::set(const KEY &key, TYPE item) {
	void *pval = NULL;
	int rc;

	// check if the key exists
	JHSG(pval, judy_hs, (char *) &key, sizeof(KEY));
	if (pval == NULL) {
		// add to array
		Word_t idx = 0;
		JLFE(rc, judy_l, idx);
		if (!rc)
			return false;
		JLI(pval, judy_l, idx);			// insert into the array
		if (pval == NULL)
			return false;

		*((TYPE **) pval) = new TYPE;
		**((TYPE **) pval) = item;

		// store the key -> value association
		JHSI(pval, judy_hs, (char *) &key, sizeof(KEY));
		*(Word_t *) pval = idx;

		return true;
	}
	else {
		// replace value in the array
		Word_t idx = *(Word_t *) pval;
		JLG(pval, judy_l, idx);			// insert into the array
		if (pval == NULL)
			return false;

		**((TYPE **) pval) = item;
		return true;
	}
}

template<class KEY, class TYPE>
bool Map<KEY, TYPE>::remove(const KEY &key) {
	void *pval = NULL;
	// check if the key exists
	JHSG(pval, judy_hs, (char *) &key, sizeof(KEY));
	if (pval == NULL) {
		// removing non-existent item
		return false;
	}
	else {
		bool res = true;
		Word_t idx = *(Word_t *) pval;
		free_item(idx);
		int rc;
		JLD(rc, judy_l, idx);
		res &= rc == 1;
		JHSD(rc, judy_hs, (char *) &key, sizeof(KEY));
		res &= rc == 1;
		return res;
	}
}

template<class KEY, class TYPE>
void Map<KEY, TYPE>::remove_all() {
	// free associated memory
	void *pval = NULL;
	Word_t idx = 0;
	JLF(pval, judy_l, idx);
	for (; idx != -1 && pval != NULL; ) {
		free_item(idx);
		JLN(pval, judy_l, idx);
	}

	int val;
	// clean array
	JLFA(val, judy_l);
	// clean index hash Map
	JHSFA(val, judy_hs);
}

template<class KEY, class TYPE>
void Map<KEY, TYPE>::free_item(Word_t idx) {
	void *pval = NULL;
	JLG(pval, judy_l, idx);
	if (pval != NULL) {
		TYPE *n = *((TYPE **) pval);
		delete n;
	}
}

template<class KEY, class TYPE>
Word_t Map<KEY, TYPE>::first() const {
	void *pval = NULL;
	Word_t idx = 0;
	JLF(pval, judy_l, idx);
	return pval ? idx : INVALID_IDX;
}

template<class KEY, class TYPE>
Word_t Map<KEY, TYPE>::next(Word_t idx) const {
	void *pval = NULL;
	JLN(pval, judy_l, idx);
	return pval ? idx : INVALID_IDX;
}

template<class KEY, class TYPE>
Word_t Map<KEY, TYPE>::last() const {
	void *pval = NULL;
	Word_t idx = -1;
	JLL(pval, judy_l, idx);
	return pval ? idx : INVALID_IDX;
}

template<class KEY, class TYPE>
Word_t Map<KEY, TYPE>::prev(Word_t idx) const {
	void *pval = NULL;
	JLP(pval, judy_l, idx);
	return pval ? idx : INVALID_IDX;
}

#endif
