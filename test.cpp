#include <iostream>
#include "ADS_set.h"
#include <random>
#include <vector>

int VALUES { 20 };

void test_insert(ADS_set<int, 2>& a) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dist(1, 100);

  for(size_t i { 0 }; i < VALUES; ++i) {
    a.insertKey(dist(gen));
  }
  std::cout << '\n';

}


void test_consecutive(ADS_set<int, 2>& a) {
  std::vector<int> values {1, 2, 3, 4, 5};
  for(const auto& i : values)
    a.insertKey(i);
}

int main() {
  ADS_set<int, 2> a;
  test_consecutive(a);
  a.dump();
}



