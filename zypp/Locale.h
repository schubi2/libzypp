/*---------------------------------------------------------------------\
|                          ____ _   __ __ ___                          |
|                         |__  / \ / / . \ . \                         |
|                           / / \ V /|  _/  _/                         |
|                          / /__ | | | | | |                           |
|                         /_____||_| |_| |_|                           |
|                                                                      |
\---------------------------------------------------------------------*/
/** \file	zypp/Locale.h
 *
*/
#ifndef ZYPP_LOCALE_H
#define ZYPP_LOCALE_H

#include <iosfwd>

#include "zypp/base/PtrTypes.h"

#include "zypp/LanguageCode.h"
#include "zypp/CountryCode.h"

///////////////////////////////////////////////////////////////////
namespace zypp
{ /////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  //
  //	CLASS NAME : Locale
  //
  /** */
  class Locale
  {
    friend std::ostream & operator<<( std::ostream & str, const Locale & obj );

  public:
    /** Implementation  */
    class Impl;

  public:
    /** Default ctor */
    Locale();

    /** Ctor taking a string. */
    explicit
    Locale( const std::string & code_r );

    /** Ctor taking LanguageCode and optional CountryCode. */
    Locale( const LanguageCode & language_r,
            const CountryCode & country_r = CountryCode() );

    /** Dtor */
    ~Locale();

  public:
    /** \name Locale constants. */
    //@{
    /** No or empty code. */
    static const Locale noCode;
    //@}

  public:
    /** */
    const LanguageCode & language() const;
    /** */
    const CountryCode & country() const;

    /** Return the locale code. */
    std::string code() const;

    /** Return the name made of language and country name. */
    std::string name() const;

    /** Return a fallback locale for this locale, when giving up, returns empty Locale() */
    Locale fallback() const;

  private:
    /** Pointer to implementation */
    RW_pointer<Impl> _pimpl;
  };
  ///////////////////////////////////////////////////////////////////

  /** \relates Locale Stream output */
  inline std::ostream & operator<<( std::ostream & str, const Locale & obj )
  { return str << obj.code(); }

  /** Comparison based on string value. */
  //@{
  /** \relates Locale */
  inline bool operator==( const Locale & lhs, const Locale & rhs ) {
    return( lhs.code() == rhs.code() );
  }
  /** \relates Locale */
  inline bool operator==( const std::string & lhs, const Locale & rhs ) {
    return( lhs == rhs.code() );
  }
  /** \relates Locale */
  inline bool operator==( const Locale & lhs, const std::string & rhs ) {
    return( lhs.code() == rhs );
  }

  /** \relates Locale */
  inline bool operator!=( const Locale & lhs, const Locale & rhs ) {
    return( ! operator==( lhs, rhs ) );
  }
  /** \relates Locale */
  inline bool operator!=( const std::string & lhs, const Locale & rhs ) {
    return( ! operator==( lhs, rhs ) );
  }
  /** \relates Locale */
  inline bool operator!=( const Locale & lhs, const std::string & rhs ) {
    return( ! operator==( lhs, rhs ) );
  }
  //@}

  /////////////////////////////////////////////////////////////////
} // namespace zypp
///////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////
namespace std
{ /////////////////////////////////////////////////////////////////
  /** \relates zypp::Locale Default order for std::container based on code string value.*/
  template<>
    inline bool less<zypp::Locale>::operator()( const zypp::Locale & lhs, const zypp::Locale & rhs ) const
    { return lhs.code() < rhs.code(); }
  /////////////////////////////////////////////////////////////////
} // namespace std
///////////////////////////////////////////////////////////////////
#endif // ZYPP_LOCALE_H
