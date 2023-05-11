#ifndef PROPELLER_UTIL_PROPERTY
#define PROPELLER_UTIL_PROPERTY

#include <functional>
namespace Foam
{
namespace util
{
/**
 * @brief class for handling properties that have some constraints
 * to be valid, like for a vector to have mag()>0 or for a scalar to be positive
 * Also handles if the variable is set or not
 * 
 * @tparam type 
 */
template<class type>
class property
{
private:
    type var_;
    bool set_ = false;
    List<std::function<bool(const type&)>> checkValid_;

public:

    property() = default;

    property(std::function<bool(const type&)> checkValid)
    {
        addCondition(checkValid);
    }

    void addCondition(std::function<bool(const type&)> checkValid)
    {
        checkValid_.append(checkValid);
    }
    static bool isPositive(const type& val) 
    {
        return val>0;
    }
    static bool isNonNegative(const type& val)
    {
        return val>=0;
    }
    static std::function<bool(const type&)> isGreaterThan(const type& other)
    {
        return [&other](const type& val)
        {
            return val > other;
        };
    }
    static bool hasLenght(const type& val)
    {
        return mag(val)>0;
    }

    bool isValid() const
    {
        bool valid=true;
        forAll(checkValid_,i)
        {
            valid &= checkValid_[i](var_);
        }
        return valid;
    }
    
    bool isSet() const
    {
        return set_;
    }

    bool isReady() const
    {
        return isSet() && isValid();
    }
    
    void set(const type& val)
    {
        set_=true;
        var_ = val;
    }
    const type& get() const 
    {
        return var_;
    }

    operator type() const { return var_; } 

    friend Ostream& operator<<(Ostream& os, const property<type>& p)
    {
        if(!p.isSet())
        {
            os<<"unset";
        }
        else if(!p.isValid())
        {
            os<<"invalid("<<p.get()<<")";
        }
        else
        {
            os<<p.get();
        }
        return os;
    }
};

}
}

#endif /* PROPELLER_UTIL_PROPERTY */
