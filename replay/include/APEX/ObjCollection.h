#ifndef APEX_ObjCollection_h
#define APEX_ObjCollection_h

#include <concepts>
#include <cstddef> 
#include <vector>
#include <algorithm>
#include <type_traits>

namespace APEX
{

// Type representing unique object ID 
using object_id = std::size_t;

/// This concept requires that the callable F is callable with single argument 'const T&' and returns 'bool'  
template<typename F, typename T> concept BoolPredicate = requires(F&& f, const T& t) {
    { std::forward<F>(f).operator()(t) } -> std::same_as<bool>; 
};

template<typename T> class ObjCollection {
private: 

    //check that the class given has a field 'T::id', and furthermore that it defines an equivalence operator 
    static_assert(requires(T& t) {
        { t.id } -> std::same_as<object_id&>; 
    }, "Type 'T' is missing 'object_id id' field / ObjCollection<T> is not friended"); 

    static_assert(std::is_move_assignable_v<T>,     "Type 'T' is not move-assignable"); 
    static_assert(std::is_move_constructible_v<T>,  "Type 'T' is not move-constructable"); 
    static_assert(std::is_copy_constructible_v<T>,  "Type 'T' is not copy-constructible");

    object_id counter{0}; 
    std::vector<T> objects{};

public: 

    using iterator = std::vector<T>::iterator; 
    using const_iterator = std::vector<T>::const_iterator; 

    ObjCollection() = default; 

    iterator begin() { return objects.begin(); }
    iterator end()   { return objects.end(); }

    const_iterator begin() const { return objects.begin(); }
    const_iterator end()   const { return objects.end(); }

    const_iterator cbegin() const { return objects.cbegin(); }
    const_iterator cend()   const { return objects.cend(); }

    std::size_t size() const { return objects.size(); }

    void reserve(std::size_t _size) { objects.reserve(_size); }

    /// @brief Get object with particular ID
    /// @param id id of object to fetch
    /// @return iterator pointing to object's posision, or this->end() if no matching object was found. 
    iterator get(object_id id) {
        return std::find_if(this->begin(), this->end(), [id](const T& rhs){ return id == rhs.id; }); 
    }

    /// @brief Add a new object to the collection. A ref to the newly placed object is returned. 
    /// @param ...args Constructor arguments  
    /// @return Size of collection after new item is added 
    template<typename...Args> 
        requires std::constructible_from<T,Args...> // <- this concept means that this function won't compile unless the 'Args&&... args' list of 
                                                    //    ctor arguments are a legal way to construct type 'T'
    object_id emplace(Args&&... args) {
        
        //first, check that this object is constructible from the 'args' list
        static_assert(
            requires(Args&&...args) {
                T(std::forward<Args>(args)...); 
            }, "Unable to construct new instance of type 'T' from list of args given. make sure all types match." 
        ); 

        T obj(std::forward<Args>(args)...); 
        obj.id = counter; 
        ++counter; 

        objects.emplace_back( std::move(obj) );
        return objects.back().id; 
    }

    /// @brief Check if matching object exists in list 
    /// @param obj object to check for match to
    /// @return 'true' if matching object is found, 'false' otherwise 
    bool exists_in_list(const T& obj) const {
        auto it = std::find_if(this->begin(), this->end(), [&obj](const T& rhs){ return obj.id == rhs.id; }); 
        return (it != this->end()); 
    }

    /// @brief Removes all elements which satisfy a user-supplied predicate 
    /// @tparam F Predicate expression 
    /// @param predicate The predicate expresion. must accept 'const T&' as argument, and return bool. removes all elements which evaluate 'predicate(elem)==true'. 
    /// @return number of elements removed. 
    template<typename F> 
        requires BoolPredicate<F,T>  // <- see comment on def. of this concept above to see what conditions it imposes. 
    std::size_t remove_if(F&& predicate) {

        std::size_t size_before = this->size(); 
        
        // this shifts all elements which evaluate 'predicate()==false' to the front, and moves all elements for which 
        // 'predicate()==true' to the front. return a new iterator pointing to the end of the range of elements for which
        // 'predicate()==false' (end of all kept elements)
        auto new_end = std::remove_if(this->begin(), this->end(), std::forward<F>(predicate)); 

        objects.erase(new_end, this->end()); 

        return size_before - this->size(); 
    }


    /// @brief Erase matching object from list, if it exists. 
    /// @param obj instance to erase
    /// @return 'true' if a matching item was removed, 'false' otherwise
    bool remove_item(const T& obj) {
        std::size_t n_removed = remove_if([&obj](const T& rhs){ return obj.id == rhs.id; }); 
        return n_removed > 0; 
    }
}; 

}

#endif