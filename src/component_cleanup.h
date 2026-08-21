/*
 * Author: Marcos Jimenez
 * Modification date: 21/08/2026
 *
 * Scope guards for factory-created model components.
 */

#ifndef LATENT_COMPONENT_CLEANUP_H
#define LATENT_COMPONENT_CLEANUP_H

template <typename T>
class pointer_vector_guard {

public:

  explicit pointer_vector_guard(std::vector<T*>& pointers):
    pointers_(pointers) {}

  pointer_vector_guard(const pointer_vector_guard&) = delete;
  pointer_vector_guard& operator=(const pointer_vector_guard&) = delete;

  ~pointer_vector_guard() {

    for(T*& pointer : pointers_) {
      delete pointer;
      pointer = nullptr;
    }

  }

private:

  std::vector<T*>& pointers_;

};

template <typename T>
class nested_pointer_vector_guard {

public:

  explicit nested_pointer_vector_guard(
      std::vector<std::vector<T*>>& pointers):
    pointers_(pointers) {}

  nested_pointer_vector_guard(const nested_pointer_vector_guard&) = delete;
  nested_pointer_vector_guard& operator=(
    const nested_pointer_vector_guard&) = delete;

  ~nested_pointer_vector_guard() {

    for(std::vector<T*>& group : pointers_) {

      for(T*& pointer : group) {
        delete pointer;
        pointer = nullptr;
      }

    }

  }

private:

  std::vector<std::vector<T*>>& pointers_;

};

#endif
