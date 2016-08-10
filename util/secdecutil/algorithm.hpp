#ifndef SecDecUtil_algorithm_hpp_included
#define SecDecUtil_algorithm_hpp_included

#include "series.hpp"

namespace secdecutil {

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /*

     Applies func to each base_type element of the nest and returns a new nest containing the result
     func may modify the elements of the original nest

     */
    template<typename out_base_type, typename in_base_type, typename T>
    struct transform_impl
    {
        using base_type = T;
        using new_type = out_base_type;
        T& nest;
        const std::function<out_base_type(in_base_type)>& func;
        new_type apply()
        {
            return func(nest);
        };
        transform_impl(T& nest, const std::function<out_base_type(in_base_type)>& func): nest(nest), func(func) {};
    };

    template<typename out_base_type, typename in_base_type, template<typename... > class E, typename T>
    struct transform_impl<out_base_type,in_base_type,E<T>>
    {
        using base_type = typename transform_impl<out_base_type,in_base_type,T>::base_type;
        using new_type = E<typename transform_impl<out_base_type,in_base_type,T>::new_type>;
        E<T>& nest;
        const std::function<out_base_type(in_base_type)>& func;
        new_type apply()
        {
            new_type content;
            for ( auto& element : nest )
                content.push_back( transform_impl<out_base_type,in_base_type,T>(element,func).apply() );
            return content;
        };
        transform_impl(E<T>& nest, const std::function<out_base_type(in_base_type)>& func): nest(nest), func(func) {};
    };

    // Specialisation for secdecutil::Series
    template<typename out_base_type, typename in_base_type, typename T>
    struct transform_impl<out_base_type,in_base_type,secdecutil::Series<T>>
    {
        using base_type = typename transform_impl<out_base_type,in_base_type,T>::base_type;
        using new_type = secdecutil::Series<typename transform_impl<out_base_type,in_base_type,T>::new_type>;
        secdecutil::Series<T>& nest;
        const std::function<out_base_type(in_base_type)>& func;
        new_type apply()
        {
            std::vector<typename transform_impl<out_base_type,in_base_type,T>::new_type> content;
            for ( auto& element : nest )
                content.push_back( transform_impl<out_base_type,in_base_type,T>(element,func).apply() );
            return secdecutil::Series<typename transform_impl<out_base_type,in_base_type,T>::new_type>
            (
             nest.get_order_min(),
             nest.get_order_max(),
             content,
             nest.get_truncated_above()
             );
        };
        transform_impl(secdecutil::Series<T>& nest, const std::function<out_base_type(in_base_type)>& func): nest(nest), func(func) {};
    };

    template<typename out_base_type, typename in_base_type, typename T>
    auto transform(T& nest, const std::function<out_base_type(in_base_type)>& func)
    -> decltype( transform_impl<out_base_type,in_base_type,T>(nest, func).apply() )
    {
        return transform_impl<out_base_type,in_base_type,T>(nest, func).apply();
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /*

     Applies func to each base_type element of the nest and returns a new nest containing the result
     func may not modify the elements of the original nest

     */
    template<typename out_base_type, typename in_base_type, typename T>
    struct const_transform_impl
    {
        using base_type = T;
        using new_type = out_base_type;
        const T& nest;
        const std::function<out_base_type(in_base_type)>& func;
        new_type apply() const
        {
            return func(nest);
        };
        const_transform_impl(const T& nest, const std::function<out_base_type(in_base_type)>& func): nest(nest), func(func) {};
    };

    template<typename out_base_type, typename in_base_type, template<typename... > class E, typename T>
    struct const_transform_impl<out_base_type,in_base_type,E<T>>
    {
        using base_type = typename const_transform_impl<out_base_type,in_base_type,T>::base_type;
        using new_type = E<typename const_transform_impl<out_base_type,in_base_type,T>::new_type>;
        const E<T>& nest;
        const std::function<out_base_type(in_base_type)>& func;
        new_type apply() const
        {
            new_type content;
            for ( const auto& element : nest )
                content.push_back( const_transform_impl<out_base_type,in_base_type,T>(element,func).apply() );
            return content;
        };
        const_transform_impl(const E<T>& nest, const std::function<out_base_type(in_base_type)>& func): nest(nest), func(func) {};
    };

    // Specialisation for secdecutil::Series
    template<typename out_base_type, typename in_base_type, typename T>
    struct const_transform_impl<out_base_type,in_base_type,secdecutil::Series<T>>
    {
        using base_type = typename const_transform_impl<out_base_type,in_base_type,T>::base_type;
        using new_type = secdecutil::Series<typename const_transform_impl<out_base_type,in_base_type,T>::new_type>;
        const secdecutil::Series<T>& nest;
        const std::function<out_base_type(in_base_type)>& func;
        new_type apply() const
        {
            std::vector<typename const_transform_impl<out_base_type,in_base_type,T>::new_type> content;
            for ( const auto& element : nest )
                content.push_back( const_transform_impl<out_base_type,in_base_type,T>(element,func).apply() );
            return secdecutil::Series<typename const_transform_impl<out_base_type,in_base_type,T>::new_type>
            (
             nest.get_order_min(),
             nest.get_order_max(),
             content,
             nest.get_truncated_above()
             );
        };
        const_transform_impl(const secdecutil::Series<T>& nest, const std::function<out_base_type(in_base_type)>& func): nest(nest), func(func) {};
    };

    template<typename out_base_type, typename in_base_type, typename T>
    auto transform(const T& nest, const std::function<out_base_type(in_base_type)>& func)
    -> decltype( const_transform_impl<out_base_type,in_base_type,T>(nest, func).apply() )
    {
        return const_transform_impl<out_base_type,in_base_type,T>(nest, func).apply();
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /*

     Applies func to each base_type element of the nest
     func may modify the elements of nest

     */
    template<typename in_base_type, typename T>
    struct void_transform_impl
    {
        T& nest;
        const std::function<void(in_base_type)>& func;
        void apply()
        {
            func(nest);
        };
        void_transform_impl(T& nest, const std::function<void(in_base_type)>& func): nest(nest), func(func) {};
    };

    template<typename in_base_type, template<typename... > class E, typename T>
    struct void_transform_impl<in_base_type,E<T>>
    {
        E<T>& nest;
        const std::function<void(in_base_type)>& func;
        void apply()
        {
            for ( auto& element : nest )
                void_transform_impl<in_base_type,T>(element,func).apply();
        };
        void_transform_impl(E<T>& nest, const std::function<void(in_base_type)>& func): nest(nest), func(func) {};
    };

    template<typename in_base_type, typename T>
    void transform(T& nest, const std::function<void(in_base_type)>& func)
    {
        void_transform_impl<in_base_type,T>(nest, func).apply();
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /*

     Applies func to each base_type element of the nest
     func may not modify the elements of nest (side-effects only)

     */
    template<typename in_base_type, typename T>
    struct void_const_transform_impl
    {
        const T& nest;
        const std::function<void(in_base_type)>& func;
        void apply() const
        {
            func(nest);
        };
        void_const_transform_impl(const T& nest, const std::function<void(in_base_type)>& func): nest(nest), func(func) {};
    };

    template<typename in_base_type, template<typename... > class E, typename T>
    struct void_const_transform_impl<in_base_type,E<T>>
    {
        const E<T>& nest;
        const std::function<void(in_base_type)>& func;
        void apply() const
        {
            for ( auto& element : nest )
                void_const_transform_impl<in_base_type,T>(element,func).apply();
        };
        void_const_transform_impl(const E<T>& nest, const std::function<void(in_base_type)>& func): nest(nest), func(func) {};
    };

    template<typename in_base_type, typename T>
    void transform(const T& nest, const std::function<void(in_base_type)>& func)
    {
        void_const_transform_impl<in_base_type,T>(nest, func).apply();
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////

}

#endif
