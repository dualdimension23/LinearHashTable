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
  size_t _allocatedSize { 0 };
  // For linear hashing, one needs to keep track of the round number
  size_t _roundNumber { 0 };
  // As well as "next to split", it points to the bucket, which will be split next
  size_t _nextToSplit { 0 };
  // Further, we need a parameter which keeps track of the amount of keys stored in the set
  size_t _elements { 0 };



public:
  ADS_set() {};
  ADS_set(std::initializer_list<key_type> ilist) : ADS_set() {
    insert(ilist);
  }
  template<typename InputIt> ADS_set(InputIt first, InputIt last) : ADS_set() {
    insert(first, last);
  }
  ADS_set(const ADS_set &other) {
    operator=(other);
  }
  ~ADS_set() {
    for(size_t index { 0 }; index < _size; ++index)
      delete _buckets[index];
    delete[] _buckets;
  }

  ADS_set &operator=(const ADS_set &other) {
    if (this == &other)
      return *this;
    clear();
    for (const auto& i : other) {
      insertKey(i);
    }
    return *this;
  }

  ADS_set &operator=(std::initializer_list<key_type> ilist) {
    clear();
    if(!ilist.size()) return *this;
    insert(ilist);
    return *this;
  }

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

  iterator find(const key_type &key) const {
    if (!_elements || !count(key)) return end();
    size_t index { getIndex(key) };
    int foundBucketIndex { -1 };
    Bucket* bucket { _buckets[index] };
    bool isPrimary { true };
    bool breakWhile { false };
    do {
      if(!isPrimary) bucket = bucket->_overflow;
      isPrimary = false;
      for(size_t bucketIndex { 0 }; bucketIndex < bucket->_size; ++bucketIndex) {
        if(key_equal{}(bucket->_records[bucketIndex], key)) {
          foundBucketIndex = bucketIndex;
          breakWhile = true;
          break;
        }
      }
      if(breakWhile) break;
    } while(bucket->_overflow);
    return Iterator(this, bucket, &key, index, foundBucketIndex);
  }

  void clear() {
    clearTable();
    setClassVariablesToDefault();
  }

  void clearTable() {
    if(!_buckets) return;
    for(size_t index { 0 }; index < _size; ++index)
      delete _buckets[index];
    delete[] _buckets;
  }

  void setClassVariablesToDefault() {
    _buckets = nullptr;
    _size = _allocatedSize = _roundNumber = _nextToSplit = _elements = 0;
  }

  void swap(ADS_set &other) {
    // Swap Idiom
    std::swap(_buckets, other._buckets);
    std::swap(_nextToSplit, other._nextToSplit);
    std::swap(_elements, other._elements);
    std::swap(_allocatedSize, other._allocatedSize);
    std::swap(_roundNumber, other._roundNumber);
    std::swap(_size, other._size);
  }

  void insert(std::initializer_list<key_type> ilist) {
    insert(std::begin(ilist), std::end(ilist));
  }
  std::pair<iterator,bool> insert(const key_type &key) {
    if(!_size) initializeSet();
    if(count(key)) return std::make_pair(find(key), false);
    insertKey(key);
    return std::make_pair(Iterator(this, findOverflowBucket(key), &key, findIndex(key), findBucketIndex(key)), true);
  }

  Bucket* findOverflowBucket(const key_type& key) const {
    if(!_elements) return nullptr;
    size_t index { getIndex(key) };
    Bucket* bucket { _buckets[index] };
    bool isPrimary { true };
    do {
      if(!isPrimary) bucket = bucket->_overflow;
      isPrimary = false;
      for(size_t bucketIndex { 0 }; bucketIndex < bucket->_size; ++bucketIndex) {
        if(key_equal{}(bucket->_records[bucketIndex], key))
          return bucket;
      }
    } while(bucket->_overflow);
    return nullptr;
  }

  int findIndex(const key_type& key) const {
    if(!_elements) return -1;
    size_t index { getIndex(key) };
    Bucket* bucket { _buckets[index] };
    bool isPrimary { true };
    do {
      if(!isPrimary) bucket = bucket->_overflow;
      isPrimary = false;
      for(size_t bucketIndex { 0 }; bucketIndex < bucket->_size; ++bucketIndex) {
        if(key_equal{}(bucket->_records[bucketIndex], key))
          return index;
      }
    } while(bucket->_overflow);
    return -1;
  }

  int findBucketIndex(const key_type& key) const {
    if(!_elements) return -1;
    size_t index { getIndex(key) };
    Bucket* bucket { _buckets[index] };
    bool isPrimary { true };
    do {
      if(!isPrimary) bucket = bucket->_overflow;
      isPrimary = false;
      for(size_t bucketIndex { 0 }; bucketIndex < bucket->_size; ++bucketIndex) {
        if(key_equal{}(bucket->_records[bucketIndex], key))
          return bucketIndex;
      }
    } while(bucket->_overflow);
    return -1;
  }

  template<typename InputIt> void insert(InputIt first, InputIt last) {
    for (auto it { first }; it != last; ++it)
      insertKey(*it);
  }

  void insertKey(const key_type& key, bool needsSplit = true) {
    // For insertion there are serveral cases to be handled.
    // First we need to check if there is even memory allocated for the set
    if (!_size) allocateMemory();
    // Second, if the key is already in the set, we do not have to insert it again
    if (count(key)) return;
    
    // We need to find the bucket where the key shall be inserted
    size_t index { getIndex(key) };
    Bucket* bucket { _buckets[index] };

    // Bucket is now found, for safety we check if the bucket still has overflow buckets,
    // since we want to insert the key at the very end
    while (bucket->_overflow) bucket = bucket->_overflow;

    // We are now at the very end
    // ( 1 ) If the bucket is full, we have to create a new overflow bucket and rehash the table
    // ( 2 ) Else if the bucket is not full, we just insert the key

    // ( 1 ) 
    if (bucket->isFull()) {
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
      bucket->_records[bucket->_size] = key;
      ++bucket->_size;
      if(needsSplit)
        ++_elements;
    }
  }

  void split() {
    // An overflow occured, we want to increase the domain and rehash the keys of the bucket where nextToSpit is pointing at
    // There are 2 cases to handle:
    // ( 1 ) We still have enough allocatedSize to just increase the shown size (_size) 
    if (_size + 1 <= _allocatedSize) enoughSize();
    // ( 2 ) If we just increase _size we would get out of bounds => we have to create a new table with more allocatedSize
    else
      allocateMemory();
  }

  void redistributeBucket(Bucket bucket) {
    // Redistribution of keys after split has happened
    // First, we need to increment nextToSplit, otherwise keys would be inserted into the bucket they are currently stored in
    ++_nextToSplit;
    bool isPrimary { true };
    do {
      if(!isPrimary) {
        bucket = *bucket._overflow;
      }
      isPrimary = false;
      for(size_t index { 0 }; index < bucket._size; ++index) 
        insertKey(bucket._records[index], false);
    } while(bucket._overflow);
  }


  void enoughSize() {
    // In split case ( 1 ) was the case
    // The next bucket is created, after that we need to increment _size
    _buckets[_size++] = new Bucket;

    // Access our splitBucket
    Bucket* splitBucket { _buckets[_nextToSplit] };
    // Store values in temporaryBucket
    Bucket temporaryBucket { *splitBucket };
    Bucket emptyBucket {};
    *splitBucket = emptyBucket;
    // The values of this splitBucket still need to be redistributed on the now (1 larger domain) table
    redistributeBucket(temporaryBucket);

    // Finally, we have to check if the round is done
    checkRound();
  }

  void initializeSet() {
    // If the table has no memory allocated, we need to increase _allocatedSize & _size
    // further we create the first bucket
    _allocatedSize = ++_size;
    _buckets = new Bucket*[_allocatedSize];
    _buckets[0] = new Bucket;
  }


  void allocateMemory() {
    // In split case ( 2 ) was the case
    // The table is simply too small, we can not increase the domain so easily
    // Lets save the current size into a variable
    size_t oldSize { _size };
    // We have to create a new set with a even larger domain ( lets say 2 times as large ) 
    size_t newSize { (2 * _allocatedSize) + 1 };
    // First lets put our old set into a temporary set
    Bucket** temporaryBuckets { _buckets };
    Bucket** newBuckets { new Bucket*[newSize] };

    // Okay, that's done. We can now set out actual set (_buckets) equal to the new set (newBuckets)
    _buckets = newBuckets;
    _allocatedSize = newSize;

    // We now have 
    //    ( 1 ) An empty set (_buckets) twice as big as our old version
    //    ( 2 ) and a copy (temporaryBuckets) of our old version

    // Remember, allocateMemory is initialized in split()
    // so we have to increase the shown size (_size) by one
    ++_size;
    if(!oldSize) {
      _buckets[0] = new Bucket;
      return;
    }

    // Now we just have to copy the buckets from temporary into our fresh new table (_buckets)
    // But ATTENTION: the splitBucket (bucket where nextToSplit points to) has to be redistributed

    // Now we can copy every bucket BUT the splitBucket into our new set
    for(size_t index { 0 }; index < oldSize; ++index) {
      if (index != _nextToSplit) {
        // We access both buckets
        // the depricated (old) bucket
        // and the empty (new) bucket
        // We literally just copy over the bucket
        newBuckets[index] = temporaryBuckets[index];
      }
    }
    
    _buckets[_nextToSplit] = new Bucket;
    _buckets[oldSize] = new Bucket;

    // We have to access our splitBucket again
    Bucket* splitBucket { temporaryBuckets[_nextToSplit] };
    redistributeBucket(*splitBucket);



    // Okay, so far so good. Time for clean up
    delete temporaryBuckets[_nextToSplit - 1];
    delete[] temporaryBuckets;
    
    // As for the "enoughSize" case, we need to check if the round is done
    checkRound();
  }


  void checkRound() {
    // We come right after a splitting of the set has happened
    // Lets check if the round is over
    if (_nextToSplit == factor(_roundNumber)) {
      // If so, we need to reset "nextToSplit" and enter the next round 
      _nextToSplit = 0;
      ++_roundNumber;
      }
  }


  size_type erase(const key_type &key) {
    // If the key is not in the set, we return 0
    if(!count(key)) return 0;
    // Get right index and access bucket
    size_t index { getIndex(key) };
    Bucket* bucket { _buckets[index] };
    bool isPrimary { true };

    do {
      if(!isPrimary) bucket = bucket->_overflow;
      isPrimary = false;
      for(size_t index { 0 }; index < bucket->_size; ++index) {
        if(key_equal{}(bucket->_records[index], key)) {
          if(index == bucket->_size) {
          } else {
            for(size_t moveIndex { index + 1 }; moveIndex < bucket->_size; ++moveIndex)
              bucket->_records[moveIndex - 1] = bucket->_records[moveIndex];
          }
          --_elements;
          --bucket->_size;
          return 1;
        }
      }
    } while(bucket->_overflow);
    return 1;
  }

  const_iterator begin() const {
    if(!_elements || !_size) return end();
    size_t index { 0 };
    Bucket* bucket { _buckets[index] };
    while(!bucket->_size && index < _size) {
      if(!bucket->_overflow)
        bucket = _buckets[++index];
      else
        bucket = bucket->_overflow;
    }
    Iterator it = Iterator(this, bucket, bucket->_records, index, 0);
    return it;
  }
  const_iterator end() const {
    return Iterator(this);
  }

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
      o << index << " Size: [" << bucket->_size << "]: |";
      bool isPrimary { true };
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

  friend bool operator==(const ADS_set &lhs, const ADS_set &rhs) {
    if(lhs.size() != rhs.size()) return false;
    for(size_t index { 0 }; index < lhs._size; ++index) {
      Bucket* bucket { lhs._buckets[index] };
      bool isPrimary { true };
      do {
        if(!isPrimary) bucket = bucket->_overflow;
        isPrimary = false;
        for(size_t bucketIndex { 0 }; bucketIndex < bucket->_size; ++bucketIndex) {
          if(!rhs.count(bucket->_records[bucketIndex])) return false;
      }
      } while(bucket->_overflow);
    }
    return true;
  }

  friend bool operator!=(const ADS_set &lhs, const ADS_set &rhs) {
    return !(lhs == rhs);
  }
};

template <typename Key, size_t N>
class ADS_set<Key,N>::Iterator {
public:
  using value_type = Key;
  using difference_type = std::ptrdiff_t;
  using reference = const value_type &;
  using pointer = const value_type *;
  using iterator_category = std::forward_iterator_tag;

private:
  const ADS_set* _set;
  const key_type* _key;
  int _index;
  int _bucketIndex;
  Bucket* _bucket;

public:
  Iterator() : _set { nullptr },
    _key { nullptr },
    _bucket { nullptr },
    _index { -1 },
    _bucketIndex { -1 }
  {}

  explicit Iterator(const ADS_set* set, Bucket* bucket, const key_type* key, int index, int bucketIndex) :
           _key { key },
           _set { set },
           _bucket { bucket },
           _index { index },
           _bucketIndex { bucketIndex }
  {}

  explicit Iterator(const ADS_set* set) :
           _key { nullptr },
           _set { set },
           _bucket { nullptr },
           _index { -1 },
           _bucketIndex { -1 }
  {}

  void outOfBounds() {
    _key = nullptr;
    _bucket = nullptr;
    _index = -1;
    _bucketIndex = -1;
  }
  reference operator*() const { return *_key; }
  pointer operator->() const { return _key; }
  Iterator &operator++() {
    // First, we need to check if the setpointer is a nullptr
    if(!_set) return *this;
    ++_bucketIndex;
    while(_bucket->_size == static_cast<size_t> (_bucketIndex)){
      if(!_bucket->_overflow) {
        if(static_cast<size_t>(_index + 1) == _set->_size) {
          outOfBounds();
          return *this;
        }
        ++_index;
        _bucket = _set->_buckets[_index];
      } else {
        _bucket = _bucket->_overflow;
      }
      _bucketIndex = 0;
    }
    _key = _bucket->_records + _bucketIndex;
    return *this;
  }

  Iterator operator++(int) {
    iterator temporaryIterator { *this };
    operator++();
    return temporaryIterator;
  }
  friend bool operator==(const Iterator &lhs, const Iterator &rhs) {
    return lhs._set == rhs._set && lhs._bucket == rhs._bucket && lhs._index == rhs._index && lhs._bucketIndex == rhs._bucketIndex;
  }

  friend bool operator!=(const Iterator &lhs, const Iterator &rhs) {
    return !(lhs == rhs);
  }
};

template <typename Key, size_t N> void swap(ADS_set<Key,N> &lhs, ADS_set<Key,N> &rhs) { lhs.swap(rhs); }

#endif // ADS_SET_H
