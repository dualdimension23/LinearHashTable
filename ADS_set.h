#ifndef ADS_SET_H
#define ADS_SET_H

#include <functional>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include "vector"

template <typename Key, size_t N = 2>
class ADS_set {
public:
  class Iterator;
  using value_type = Key;
  using key_type = Key;
  using reference = key_type &;
  using const_reference = const key_type &;
  using size_type = size_t;
  using difference_type = std::ptrdiff_t;
  using iterator = Iterator;
  using const_iterator = Iterator;
  using key_compare = std::less<key_type>;   // B+-Tree
  using key_equal = std::equal_to<key_type>; // Hashing
  using hasher = std::hash<key_type>;        // Hashing

private:
  // Bucket datastructure
  struct Bucket {
    // Variables
    key_type _records[N];
    Bucket* _overflow { nullptr };
    size_t _size { 0 };
    bool _isPrimary { true };
    // Constructors
    Bucket() {}
    Bucket(bool isPrimary): _isPrimary { isPrimary } {}
    // Methods
    bool isFull() { return _size == N; };
    ~Bucket() {
      if(_overflow) delete _overflow;
    }
  };

  // The set consists of an array of pointers to bucket pointers
  Bucket** _buckets { nullptr };
  // It has a "shown size" - the number of buckets exposed to the functionality
  size_t _size { 0 };
  // As well as an "actual size" - the size which is allocated
  size_t _allocatedSize { 5 };
  // For linear hashing, one needs to keep track of the round number
  size_t _roundNumber { 0 };
  // As well as "next to split", it points to the bucket, which will be split next
  size_t _nextToSplit { 0 };
  // Further, we need a parameter which keeps track of the amount of keys stored in the set
  size_t _elements { 0 };



public:
  ADS_set() {};
  //ADS_set(std::initializer_list<key_type> ilist);
  //template<typename InputIt> ADS_set(InputIt first, InputIt last);
  //ADS_set(const ADS_set &other);
  ~ADS_set() {
    for(size_t index { 0 }; index < _allocatedSize; ++index)
      delete _buckets[index];
    delete[] _buckets;
  }

  //ADS_set &operator=(const ADS_set &other);
  //ADS_set &operator=(std::initializer_list<key_type> ilist);

  size_t factor(size_t roundNumber) const {
      size_t base { 2 };
      size_t exponent { roundNumber };
      size_t result { 1 };
      for (size_t i { 0 }; i < exponent; ++i) result *= base;
      return result;
  }

  size_t firstHashFunction(const key_type& key) const {
    // Hash function used if the bucket was not split in this round
    return hasher{}(key) % factor(this->_roundNumber);
  }

  size_t secondHashFunction(const key_type& key) const {
    // Hash function used if the bucket was split in this round
    return hasher{}(key) % factor(this->_roundNumber + 1);
  }
  
  size_t getIndex(const key_type& key) const {
    // This function will be used throughout the whole program
    // it returns the index of the bucket where a specific key is thought to be
    size_t index { firstHashFunction(key) };
    if (index < _nextToSplit)
      index = secondHashFunction(key);
    return index;
  }

  size_type size() const {
    // Returns the current amount of elements, stored in the set
    return this->_elements;
  }

  bool empty() const {
    // Returns true if the set is empty, else false
    return (!this->_elements);
  }

  Bucket* findKey(const key_type& key) const {
    // Returns the bucket where the key is stored, or nullptr if the key is not found
    if (!this->_size) return nullptr;
    size_t index = { getIndex(key) };
    Bucket* bucket { _buckets[index] };
    bool isPrimary = true;

    do {
      if (!isPrimary) bucket = bucket->_overflow;
      isPrimary = false;
      for(size_t index { 0 }; index < bucket->_size; ++index) {
        if (key_equal{}(bucket->_records[index], key)) return bucket; 
      }
    } while (bucket->_overflow);
    return nullptr;
  }

  size_type count(const key_type &key) const {
    // Returns 1 if the key is found, else 0
    return !!findKey(key);
  }
  //iterator find(const key_type &key) const;

  //void clear();
  //void swap(ADS_set &other);

  void insert(std::initializer_list<key_type> ilist) {
    insert(std::begin(ilist), std::end(ilist));
  }
  //std::pair<iterator,bool> insert(const key_type &key);
  template<typename InputIt> void insert(InputIt first, InputIt last) {
    for (auto it { first }; it != last; ++it)
      insertKey(*it);
  }
  void insertKey(const key_type& key, bool needsSplit = true) {
    // For insertion there are serveral cases to be handled.
    // First we need to check if there is even memory allocated for the set
    // Second, if the key is already in the set, we do not have to insert it again
    if (!_size) allocateMemory();
    if (count(key)) return;
    
    // We need to find the bucket where the key shall be inserted
    size_t index { getIndex(key) };
    std::cout << "Inserting " << key << " at index " << index << '\n';
    Bucket* bucket { _buckets[index] };
    std::cout << "Bucket found...\n";

    // Bucket is now found, for safety we check if the bucket still has overflow buckets,
    // since we want to insert the key at the very end
    while (bucket->_overflow) bucket = bucket->_overflow;

    // We are now at the very end
    // ( 1 ) If the bucket is full, we have to create a new overflow bucket and rehash the table
    // ( 2 ) Else if the bucket is not full, we just insert the key

    // ( 1 ) 
    if (bucket->isFull()) {
      std::cout << "Bucket is full...\n";
      // We create a new overflow bucket, with isPrimary set to false
      bucket->_overflow = new Bucket(false);
      // Go into the overflow bucket
      bucket = bucket->_overflow;
      // Insert value
      bucket->_records[0] = key;
      // Set properties
      ++bucket->_size;
      if(needsSplit) {
        ++_elements;
        split();
      }
    } else {
      // ( 2 ) Okay, if the bucket is not full, we just need to find the right spot to insert
      // which is: [bucket->_size]
      std::cout << "Bucket is not full...\n";
      bucket->_records[bucket->_size] = key;
      std::cout << "Set key at index " << bucket->_size << " to " << key << '\n';
      ++bucket->_size;
      if(needsSplit)
        ++_elements;
    }
  }

  void split() {
    // An overflow occured, we want to increase the domain and rehash the keys of the bucket where nextToSpit is pointing at
    // There are 2 cases to handle:
    // ( 1 ) We still have enough allocatedSize to just increase the shown size (_size) 
    // ( 2 ) If we just increase _size we would get out of bounds => we have to create a new table with more allocatedSize
    //
    // ( 1 ) 
    if (_size + 1 < _allocatedSize) enoughSize();
    // ( 2 )
    else
      allocateMemory();
  }

  void redistributeBucket(Bucket bucket) {
    // Redistribution of keys after split has happened
    // First, we need to increment nextToSplit, otherwise keys would be inserted into the bucket they are currently stored in
    std::cout << "Redistributing... \n";
    dump();
    ++_nextToSplit;
    bool isPrimary { true };
    std::cout << "Temporary Bucket Size: " << bucket._size << '\n';
    do {
      if(!isPrimary) {
        bucket = *bucket._overflow;
        std::cout << "Another overflow...\n";
      }
      isPrimary = false;
      for(size_t index { 0 }; index < bucket._size; ++index) 
        insertKey(bucket._records[index], false);
    } while(bucket._overflow);
  }


  void enoughSize() {
    // In split case ( 1 ) was the case
    std::cout << "Entered enoughSize...\n_nextToSplit: " << _nextToSplit << '\n';
    // Lets get the bucket, which has to be rehashed
    Bucket* splitBucket { _buckets[_nextToSplit] };
    Bucket temporaryBucket { *splitBucket };
    Bucket emptyBucket {};
    *splitBucket = emptyBucket;
    // The values of this splitBucket still need to be redistributed on the now (1 larger domain) table
    std::cout << "Redistributing bucket " << _nextToSplit << "...\n";
    redistributeBucket(temporaryBucket);

    // We are now done redistributing our keys
    // Increase domain
    ++_size;
    // Finally, we have to check if the round is done
    checkRound();
  }

  void allocateMemory() {
    // In split case ( 2 ) was the case
    // The table is simply too small, we can not increase the domain so easily
    // Lets save the current size into a variable
    size_t old_size { _size };
    // We have to create a new set with a even larger domain ( lets say 2 times as large ) 
    size_t new_size { (2 * _allocatedSize) + 1 };
    // First lets put our old set into a temporary set
    std::cout << "Allocating memory...\n";
    Bucket** temporaryBuckets { _buckets };
    Bucket** newBuckets { new Bucket*[new_size] };
    for(size_t index { 0 }; index < new_size; ++index) newBuckets[index] = new Bucket;
    std::cout << "Created new table...\n";

    // Okay, that's done. We can now set out actual set (_buckets) equal to the new set (newBuckets)
    _buckets = newBuckets;
    _allocatedSize = new_size;
    
    // We now have 
    //    ( 1 ) An empty set (_buckets) twice as big as our old version
    //    ( 2 ) and a copy (temporaryBuckets) of our old version

    // Remember, allocateMemory is initialized in split(), so we have to increase the shown size (_size) by one
    // - thats the whole point of it
    ++_size;
    if(_size == 1) {
      _buckets[0] = new Bucket;
      return;
    }

    // Now we just have to copy the buckets from temporary into our fresh new table (_buckets)
    // But ATTENTION: the splitBucket (bucket where nextToSplit points to) has to be redistributed

    // Now we can copy every bucket BUT the splitBucket into our new set
    for(size_t index { 0 }; index < old_size; ++index) {
      if (index != _nextToSplit) {
        // We access both buckets
        // the depricated (old) bucket
        // and the empty (new) bucket
        Bucket* depricatedBucket { temporaryBuckets[index] };
        Bucket* newBucket { newBuckets[index] };
        // We literally just copy over the bucket
        *newBucket = *depricatedBucket;
      }
    }
    std::cout << "Copied buckets...\n";
    
    // We have to access our splitBucket again
    Bucket* splitBucket { temporaryBuckets[_nextToSplit] };
    std::cout << "Found splitBucket...\n";

    redistributeBucket(splitBucket);
    std::cout << "Redistribution finished...\n";


    // As for the "enoughSize" case, we need to check if the round is done
    checkRound();

    // Okay, so far so good. Now we can delete the temporary table again
    for(size_t index { 0 }; index < old_size; ++index)
      delete temporaryBuckets[index];
    delete[] temporaryBuckets;
    std::cout << "Deleted temporary table...\n";
  }


  void checkRound() {
    // We come right after a splitting of the set has happened
    // Lets check if the round is over
    if (_nextToSplit == factor(_roundNumber)) {
      // If so, we need to reset "nextToSplit" and enter the next round 
      _nextToSplit = 0;
      ++_roundNumber;
      std::cout << "Round " << _roundNumber << '\n';
      }
  }


  //size_type erase(const key_type &key);

  //const_iterator begin() const;
  //const_iterator end() const;

  void dump(std::ostream &o = std::cerr) const {
    // Visual Dump of the set
    // Properties
    o << "---------------------------------------------------\n";
    o << "|    _size: " << _size << "                                     \n";
    o << "|    _allocatedSize: " << _allocatedSize << "               \n";
    o << "|    _elements: " << _elements << "                         \n";
    o << "|    _roundNumber: " << _roundNumber << "                   \n";
    o << "|    _nextToSplit: " << _nextToSplit << "                   \n";
    o << "---------------------------------------------------\n\n";
    
    for(size_t index { 0 }; index < _size; ++index) {
      Bucket* bucket { _buckets[index] };
      o << index << ": |";
      bool isPrimary = true;
      do {
        if(!isPrimary) {
          bucket = bucket->_overflow;
          o << " --> "; 
        }
        isPrimary = false;
        for(size_t bucketIndex { 0 }; bucketIndex < bucket->_size; ++bucketIndex) {
          o << ' ' << bucket->_records[bucketIndex] << " |"; 
        }
      } while(bucket->_overflow);
      o << '\n';

    }
    o << "\n\n";
    
  }

  //friend bool operator==(const ADS_set &lhs, const ADS_set &rhs);
  //friend bool operator!=(const ADS_set &lhs, const ADS_set &rhs);
};

//template <typename Key, size_t N>
//class ADS_set<Key,N>::Iterator {
//public:
//  using value_type = Key;
//  using difference_type = std::ptrdiff_t;
//  using reference = const value_type &;
//  using pointer = const value_type *;
//  using iterator_category = std::forward_iterator_tag;
//
//  explicit Iterator(/* implementation-dependent */);
//  reference operator*() const;
//  pointer operator->() const;
//  Iterator &operator++();
//  Iterator operator++(int);
//  friend bool operator==(const Iterator &lhs, const Iterator &rhs);
//  friend bool operator!=(const Iterator &lhs, const Iterator &rhs);
//};

//template <typename Key, size_t N> void swap(ADS_set<Key,N> &lhs, ADS_set<Key,N> &rhs) { lhs.swap(rhs); }

#endif // ADS_SET_H
