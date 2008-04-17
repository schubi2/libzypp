#include <stdio.h>
#include <iostream>
#include <iterator>
#include <boost/test/auto_unit_test.hpp>
#include <list>

#include "zypp/ZYppFactory.h"
#include "zypp/PoolQuery.h"
#include "zypp/PoolQueryUtil.tcc"
#include "zypp/TmpPath.h"
#include "zypp/Locks.h"

#define BOOST_TEST_MODULE Locks

using std::cout;
using std::endl;
using std::string;
using namespace zypp;
using namespace boost::unit_test;

bool isLocked( const sat::Solvable & solvable )
{
  zypp::PoolItem pi( zypp::ResPool::instance().find( solvable ) );
  if( pi.status().isLocked() )
    return true;
  return false;
}

static void init_pool()
{
  Pathname dir(TESTS_SRC_DIR);
  dir += "/zypp/data/PoolQuery";

  ZYpp::Ptr z = getZYpp();
  ZConfig::instance().setSystemArchitecture(Arch("i586"));

  RepoInfo i1; i1.setAlias("factory");
  sat::Pool::instance().addRepoSolv(dir / "factory.solv", i1);
  RepoInfo i2; i2.setAlias("@System");
  sat::Pool::instance().addRepoSolv(dir / "@System.solv", i2);
}

BOOST_AUTO_TEST_CASE(pool_query_init)
{
  init_pool();
}

/////////////////////////////////////////////////////////////////////////////
//  0xx basic queries
/////////////////////////////////////////////////////////////////////////////

// default query + one search string
// q.addString("foo");
// result: all resolvables having at least one attribute matching foo
BOOST_AUTO_TEST_CASE(locks_1)
{
  cout << "****001****"  << endl;
  PoolQuery q;
  q.addString("zypper");
  Locks::instance().addLock(q);
  for_(it,q.begin(),q.end())
  {
    BOOST_CHECK(isLocked(*it));
  }
  Locks::instance().removeLock(q); //clear before next test
}

BOOST_AUTO_TEST_CASE(locks_save_load)
{
  cout << "****save/load****"  << endl;
  Pathname src(TESTS_SRC_DIR);
    src += "zypp/data/Locks/locks";
  Locks::instance().readAndApply(src);
  PoolQuery q;
  q.addString("zypper");
  for_(it,q.begin(),q.end())
  {
    BOOST_CHECK(isLocked(*it));
  }
#if 1 
  filesystem::TmpFile testfile;
  //Pathname testfile(TESTS_SRC_DIR);
    //  testfile += "/zypp/data/Locks/testlocks";
  Locks::instance().removeLock(q); 
  Locks::instance().save(testfile);
  Locks::instance().readAndApply(testfile);
  //now unlocked
  for_(it,q.begin(),q.end())
  {
    BOOST_CHECK(!isLocked(*it));
  }
  BOOST_CHECK(Locks::instance().size()==0);
#endif
}

BOOST_AUTO_TEST_CASE(locks_save_without_redundancy)
{
  cout << "****save without redundancy****"  << endl;
  PoolQuery q;
  q.addString("zypper");
  Locks& locks = Locks::instance();
  locks.addLock(q);
  locks.addLock(q);
  locks.save("/dev/null");
  BOOST_CHECK( locks.size()==1 );
  locks.addLock(q);
  locks.save("/dev/null");
  BOOST_CHECK( locks.size()==1 );
  locks.removeLock(q);
  locks.save("/dev/null");
  BOOST_CHECK( locks.size() == 0 );
}

BOOST_AUTO_TEST_CASE( locks_empty )
{
  cout << "****test and clear empty locks****"  << endl;
  PoolQuery q;
  q.addString("foo-bar-nonexist");
  Locks& locks = Locks::instance();
  locks.addLock(q);
  locks.save( "/dev/null" ); //only need merge list
  BOOST_CHECK( locks.existEmpty() );
  locks.removeEmpty();
  BOOST_CHECK( locks.size() == 0 );
}
