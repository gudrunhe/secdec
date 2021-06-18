#ifndef SecDecUtil_integrand_container_hpp_included
#define SecDecUtil_integrand_container_hpp_included

#ifdef SECDEC_WITH_CUDA
    #include <thrust/complex.h>
#endif
#include <complex>
#include <functional>
#include <memory>
#include <sys/mman.h>
#include <atomic>
#include <vector>

namespace secdecutil {

    // this error is thrown if the sign check of the deformation (contour_deformation_polynomial.imag() <= 0) fails
    struct sign_check_error : public std::runtime_error { using std::runtime_error::runtime_error; };

    // struct that contains fields that the integrand function can write information to, this gets passed to the gpu and other threads

    struct ResultInfo
    {
        std::atomic<int> filled{0};

        enum class ReturnValue { no_error, sign_check_error_contour_deformation, sign_check_error_positive_polynomial};
        ReturnValue return_value = ReturnValue::no_error;
        int signCheckId;

        void process_errors() const{
            if(return_value == ReturnValue::sign_check_error_contour_deformation){
                throw secdecutil::sign_check_error("\"contour deformation polynomial\", signCheckId=" + std::to_string(signCheckId));
            }
            else if(return_value == ReturnValue::sign_check_error_positive_polynomial){
                throw secdecutil::sign_check_error("\"positive polynomial\", signCheckId=" + std::to_string(signCheckId));
            }
        }

        void clear_errors(){
            return_value = ReturnValue::no_error;
            filled = 0;
        }

        #ifdef SECDEC_WITH_CUDA
        __device__ __host__
        #endif
        void fill_if_empty_threadsafe(const ResultInfo& result_info_new){
          #ifdef __CUDA_ARCH__
            if(not atomicCAS(reinterpret_cast<int*>(&filled), 0, 1)){
          #else
            int tempvariable = 0;
            if(std::atomic_compare_exchange_strong(&filled, &tempvariable, 1)){
          #endif
                return_value = result_info_new.return_value;
                signCheckId = result_info_new.signCheckId;
            }
        }
    };

    template<typename T, typename Args, typename Pars=double, typename Pars_extra=Pars>
    struct IntegrandContainer {
        
        template<typename T_real> struct remove_complex { using type = T_real; };
        template<typename T_real> struct remove_complex<std::complex<T_real>> { using type = T_real; };
        #ifdef SECDEC_WITH_CUDA
            template<typename T_real> struct remove_complex<thrust::complex<T_real>> { using type = T_real; };
        #endif
        
        enum Operation { add, subtract, multiply, divide };
        int number_of_integration_variables;
        std::string display_name = "INTEGRAND";
        std::vector<std::vector<Pars>> parameters;
        std::vector<std::vector<Pars*>> parameters_ptr;
        std::vector<std::vector<Pars_extra>> extra_parameters;
        std::function<T(Args,const Pars*,ResultInfo*)> integrand_with_parameters;
        std::shared_ptr<ResultInfo> result_info;

        auto get_parameters() -> decltype(parameters_ptr){return parameters_ptr;}
        auto get_extra_parameters() -> decltype(extra_parameters){return extra_parameters;}

        IntegrandContainer(const int number_of_integration_variables, const std::function<T(Args,const Pars*,ResultInfo*)>& integrand_with_parameters, std::vector<std::vector<Pars>> parameters):
            number_of_integration_variables(number_of_integration_variables), integrand_with_parameters(integrand_with_parameters), parameters(parameters), parameters_ptr()
        {
            for(int k = 0; k < parameters.size(); k++){
                parameters_ptr.push_back(std::vector<Pars*>(number_of_integration_variables));
                for(int i = 0; i < number_of_integration_variables; i++){
                    parameters_ptr[k][i] = &this->parameters[k][i];
                }
            }
            ResultInfo* result_info_raw = (ResultInfo*)mmap(NULL, sizeof(ResultInfo), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS,-1,0);
            result_info = std::shared_ptr<ResultInfo>(result_info_raw, [](ResultInfo* result_info_raw){munmap(result_info_raw,sizeof(ResultInfo));});
        }

        IntegrandContainer(const int number_of_integration_variables, const std::function<T(Args, ResultInfo*)>& integrand):
            integrand_with_parameters([integrand](Args x, const Pars* p, ResultInfo* result_info){return integrand(x,result_info);}),
            parameters(), extra_parameters() {
                this->number_of_integration_variables = number_of_integration_variables;
                ResultInfo* result_info_raw = (ResultInfo*)mmap(NULL, sizeof(ResultInfo), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS,-1,0);
                result_info = std::shared_ptr<ResultInfo>(result_info_raw, [](ResultInfo* result_info_raw){munmap(result_info_raw,sizeof(ResultInfo));});
            }

        IntegrandContainer() :
        number_of_integration_variables(0), integrand_with_parameters([](Args x, const Pars* p, ResultInfo* result_info){return T();}),
        parameters(), extra_parameters(), result_info(std::make_shared<ResultInfo>())
        {
        }

        IntegrandContainer(const IntegrandContainer& other):
        IntegrandContainer(other.number_of_integration_variables, other.integrand_with_parameters, other.parameters){
            result_info = other.result_info;
            extra_parameters = other.extra_parameters;
            display_name = other.display_name;
        }
        
        // Constructor of complex IntegrandContainer from real IntegrandContainer
        #ifdef SECDEC_WITH_CUDA
            template<typename Tc = T, typename = typename std::enable_if<std::is_same<thrust::complex<typename remove_complex<Tc>::type>,Tc>::value>::type>
        #else
            template<typename Tc = T, typename = typename std::enable_if<std::is_same<std::complex<typename remove_complex<Tc>::type>,Tc>::value>::type>
        #endif
        IntegrandContainer(const IntegrandContainer<typename remove_complex<Tc>::type,Args,Pars,Pars_extra>& other)
        {
            number_of_integration_variables = other.number_of_integration_variables;
            display_name = other.display_name;
            parameters = other.parameters;
            parameters_ptr = other.parameters_ptr;
            extra_parameters = other.extra_parameters;
            integrand_with_parameters = [other] (Args x, const Pars* p, ResultInfo* r){ return other.integrand_with_parameters(x,p,r); }; // returns T rather than remove_complex<T>
            result_info = other.result_info;
        }
        
        T operator()(Args x) const
        {
            const Pars* dataptr = parameters.size() ? parameters[0].data() : nullptr;
            return integrand_with_parameters(x, dataptr, result_info.get());
        }

        T operator()(Args x, ResultInfo* result_info) const
        {
            const Pars* dataptr = parameters.size() ? parameters[0].data() : nullptr;
            return integrand_with_parameters(x, dataptr, result_info);
        }

        void process_errors() const{
            result_info->process_errors();
        }

        void clear_errors(){
            result_info->clear_errors();
        }

        /*
         *  Helper functions
         */
        template<int operation>
        static IntegrandContainer add_subtract_multiply_or_divide(const IntegrandContainer& ic1, const IntegrandContainer& ic2)
        {
            int number_of_integration_variables = std::max(ic1.number_of_integration_variables, ic2.number_of_integration_variables);
            std::function<T(Args,const Pars*,ResultInfo*)> integrand_with_parameters;

            if (operation == add) {
                integrand_with_parameters = [ic1, ic2] (Args x, const Pars* p, ResultInfo* r) { return ic1(x,r) + ic2(x,r); };
            } else if (operation == subtract ) {
                integrand_with_parameters = [ic1, ic2] (Args x, const Pars* p, ResultInfo* r) { return ic1(x,r) - ic2(x,r); };
            } else if (operation == multiply ) {
                integrand_with_parameters = [ic1, ic2] (Args x, const Pars* p, ResultInfo* r) { return ic1(x,r) * ic2(x,r); };
            } else if ( operation == divide ) {
                integrand_with_parameters = [ic1, ic2] (Args x, const Pars* p, ResultInfo* r) { return ic1(x,r) / ic2(x,r); };
            }

            return IntegrandContainer(number_of_integration_variables, integrand_with_parameters, {});
        }

        /*
         * unary operators
         */
        IntegrandContainer operator+() const
        {
            return *this;
        };

        IntegrandContainer operator-() const
        {
            return IntegrandContainer
            (
                number_of_integration_variables,
                [this] (Args x, const Pars* p, ResultInfo* r) { return - this->integrand_with_parameters(x,p,r); }, parameters
            );
        };

        /*
         *  Compound assignment operators
         */
        IntegrandContainer& operator-=(const IntegrandContainer& ic1)
        {
            *this = *this - ic1;
            return *this;
        };

        IntegrandContainer& operator+=(const IntegrandContainer& ic1)
        {
            *this = *this + ic1;
            return *this;
        };

        IntegrandContainer& operator*=(const IntegrandContainer& ic1)
        {
            *this = *this * ic1;
            return *this;
        };
        IntegrandContainer& operator/=(const IntegrandContainer& ic1)
        {
            *this = *this / ic1;
            return *this;
        };

        /*
         *  Binary operators
         */
        friend IntegrandContainer operator-(const IntegrandContainer& ic1, const IntegrandContainer& ic2)
        {
            return add_subtract_multiply_or_divide<subtract>(ic1,ic2);
        };

        friend IntegrandContainer operator+(const IntegrandContainer& ic1, const IntegrandContainer& ic2)
        {
            return add_subtract_multiply_or_divide<add>(ic1,ic2);
        };

        friend IntegrandContainer operator*(const IntegrandContainer& ic1, const IntegrandContainer& ic2)
        {
            return add_subtract_multiply_or_divide<multiply>(ic1,ic2);
        };
        friend IntegrandContainer operator/(const IntegrandContainer& ic1, const IntegrandContainer& ic2)
        {
            return add_subtract_multiply_or_divide<divide>(ic1,ic2);
        };
    };

    namespace complex_to_real {

        #define COMPLEX_INTEGRAND_CONTAINER_TO_REAL(complex_template, OPERATION) \
        template<typename T, typename Args, typename Pars> \
        IntegrandContainer<T, Args, Pars> OPERATION(const IntegrandContainer<complex_template<T>, Args, Pars>& other) \
        { \
            IntegrandContainer<T, Args, Pars> ic; \
            ic.number_of_integration_variables = other.number_of_integration_variables; \
            ic.display_name = other.display_name; \
            ic.parameters = other.parameters; \
            ic.parameters_ptr = other.parameters_ptr; \
            ic.extra_parameters = other.extra_parameters; \
            ic.integrand_with_parameters = [other] (Args x, const Pars* p, ResultInfo* r){ return other.integrand_with_parameters(x,p,r).OPERATION(); }; \
            ic.result_info = other.result_info; \
            return ic; \
        }
        
        COMPLEX_INTEGRAND_CONTAINER_TO_REAL(std::complex,real)
        COMPLEX_INTEGRAND_CONTAINER_TO_REAL(std::complex,imag)
        #ifdef SECDEC_WITH_CUDA
          COMPLEX_INTEGRAND_CONTAINER_TO_REAL(thrust::complex,real)
          COMPLEX_INTEGRAND_CONTAINER_TO_REAL(thrust::complex,imag)
        #endif
        
        #undef COMPLEX_INTEGRAND_CONTAINER_TO_REAL

    }
}

#endif
