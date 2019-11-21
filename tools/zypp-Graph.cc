#define INCLUDE_TESTSETUP_WITHOUT_BOOST
#include "zypp/../tests/lib/TestSetup.h"
#undef  INCLUDE_TESTSETUP_WITHOUT_BOOST

#include <boost/program_options.hpp>

#include <iostream>
#include <fstream>
#include <stack>

#include <zypp/PoolQuery.h>
#include <zypp/ResObjects.h>
#include <zypp/ui/SelectableTraits.h>

using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using namespace zypp;

#define MSG cerr

///////////////////////////////////////////////////////////////////
// getopts
static std::string appname { "OOps" };

inline std::ostream & appmsg()
{ return MSG; }

inline std::ostream & appmsg( const std::string & msg_r )
{ return MSG << msg_r << endl; }

inline void apperror( const std::string & msg_r )
{ MSG << appname << ": " << msg_r << endl; }

inline void appexit( const std::string & msg_r, int exit_r ) __attribute__((noreturn));
inline void appexit( const std::string & msg_r, int exit_r )
{ apperror( msg_r ); exit( exit_r ); }

namespace zypp
{
  namespace filesystem {
    // need an istream able Pathname for parseOpts
    std::istream & operator>>( std::istream & istr, Pathname & val )
    {
      std::string p;
      istr >> p;
      val = p;
      return istr;
    }
  } // namespace filesystem
} // namespace zypp
namespace getopts
{
  namespace po = boost::program_options;
  using po::variables_map;
  variables_map parseOpts( int argc, const char * argv[] )
  {
    variables_map ret;
    try{
      po::options_description options( "OPTIONS" );
      options.add_options()
      ( "root", po::value<Pathname>()->value_name("ROOTDIR")->default_value("/"),
	"Load repos from the system located below ROOTDIR." )
      ( "PKG_SPEC", po::value<std::vector<std::string>>(),
	"Pkgspec to process")
      ( "help,h",	"Show this help message and exit" )
      ;

      po::positional_options_description positionals;
      positionals.add( "PKG_SPEC", -1 );

      po::store(
	po::command_line_parser( argc, argv )
	.options( options )
	.positional( positionals )
	.run(), ret
      );
      po::notify( ret );

      if ( ret.count( "help" ) || ! ret.count( "PKG_SPEC" ) )
      {
	std::string appname { Pathname::basename( argv[0] ) };
	cout
	<< "Usage: " << appname << " [OPTIONS] PKG_SPEC..." << endl
	<< endl
	<< options << endl;
	exit( 0 );
      }
    }
    catch ( const std::exception & excpt ) {
      appexit( excpt.what(), 1 );
    }
    return ret;
  }
} // namespace getopts

///////////////////////////////////////////////////////////////////
// poolsetup
void loadSystem( const getopts::variables_map & opts_r )
{
  const Pathname & sysRoot { opts_r["root"].as<Pathname>() };
  if ( ! PathInfo( sysRoot ).isDir() )
    return appexit( "--root requires a directory.", 1 );

  ZConfig::instance();
  sat::Pool satpool( sat::Pool::instance() );

  if ( TestSetup::isTestcase( sysRoot ) )
  {
    appmsg( str::Format("*** Load Testcase from '%s'") % sysRoot );
    TestSetup test;
    test.loadTestcaseRepos( sysRoot );
    dumpRange( appmsg(), satpool.reposBegin(), satpool.reposEnd() ) << endl;
  }
  else if ( TestSetup::isTestSetup( sysRoot ) )
  {
    appmsg( str::Format("*** Load TestSetup from '%s'") % sysRoot );
    const char * astr = getenv( "ZYPP_TESTSUITE_FAKE_ARCH" );
    if ( !astr || !*astr )
      astr = getenv( "ZYPP_ARCH" );
    if ( !astr || !*astr )
      astr = "x86_64";
    TestSetup test( sysRoot, Arch( astr ) );
    test.loadRepos();
    dumpRange( appmsg(), satpool.reposBegin(), satpool.reposEnd() ) << endl;
  }
  else
  {
    // a system
    appmsg( str::Format("*** Load system at '%s'") % sysRoot.c_str() );
    if ( false )
    {
      appmsg() << "*** load target '" << Repository::systemRepoAlias() << "'\t" << flush;
      getZYpp()->initializeTarget( sysRoot );
      getZYpp()->target()->load();
      appmsg() << satpool.systemRepo() << endl;
    }

    if ( true )
    {
      RepoManager repoManager( sysRoot );
      RepoInfoList repos = repoManager.knownRepositories();
      for_( it, repos.begin(), repos.end() )
      {
	RepoInfo & nrepo( *it );

	if ( ! nrepo.enabled() )
	  continue;

	if ( ! repoManager.isCached( nrepo ) )
	{
	  appmsg() << str::form( "*** omit uncached repo '%s' (do 'zypper refresh')", nrepo.name().c_str() ) << endl;
	  continue;
	}

	appmsg() << str::form( "*** load repo '%s'\t", nrepo.name().c_str() ) << flush;
	try
	{
	  repoManager.loadFromCache( nrepo );
	  appmsg() << satpool.reposFind( nrepo.alias() ) << endl;
	}
	catch ( const Exception & exp )
	{
	  appmsg() << exp.asString() + "\n" + exp.historyAsString() << endl;
	  appmsg() << str::form( "*** omit broken repo '%s' (do 'zypper refresh')", nrepo.name().c_str() ) << endl;
	  continue;
	}
      }
    }
  }
}

///////////////////////////////////////////////////////////////////
// the Graph
namespace depgraph
{
  using NodeId = std::string;

  ///////////////////////////////////////////////////////////////////
  /// \brief Node is a set of Solvables.
  /// It's "dependencies" is the union of it's Solvables dependencies.
  struct Node
  {
    Node( NodeId id_r )
    : _id { std::move(id_r) }
    {}

    const NodeId & id() const
    { return _id; }

  private:
    NodeId _id;
  };
  ///////////////////////////////////////////////////////////////////
  /// \brief Vertex is a dependency between 2 Nodes
  ///
  struct Vertex
  {
    Vertex( Node from_r, Node to_r )
    : _from { std::move(from_r) }
    , _to { std::move(to_r) }
    {}

  private:
    Node _from;
    Node _to;
  };
} // namespace depgraph
namespace std
{
  template <> struct hash<depgraph::Node>
  {
    size_t operator()( const depgraph::Node & el ) const
    { return std::hash<depgraph::NodeId>( el.id() ); }
  };
  template <> struct hash<depgraph::Vertex>
  {
    size_t operator()( const depgraph::Vertex & el ) const
    {

    }
  };
} //namespace std
namespace depgraph
{
  ///////////////////////////////////////////////////////////////////
  struct Graph
  {
  public:
    Graph & addNode( sat::Solvable solv_r )
    {
      _nodes.insert( Node( nodeId( solv_r ) ) );
      return *this;
    }

  public:
    /** Maps Solvable to a Node */
    NodeId nodeId( sat::Solvable solv_r ) const
    { return solv_r.ident().asString();
      return solv_r.asString(); }

  private:
    std::unordered_set<Node>   _nodes;
    std::unordered_set<Vertex> _vertices;
  };

} // namespace depgraph

///////////////////////////////////////////////////////////////////
// main
int main( int argc, const char * argv[] )
{
  appname = Pathname::basename( argv[0] );
  getopts::variables_map vm { getopts::parseOpts( argc, argv ) };
  ///////////////////////////////////////////////////////////////////
  loadSystem( vm );
  ///////////////////////////////////////////////////////////////////

  using depgraph::Graph;
  Graph g;
  for( const std::string pkgspec : vm["PKG_SPEC"].as<std::vector<std::string>>() )
  {
    PoolQuery q;
    q.setMatchGlob();
    q.setCaseSensitive( false );
    q.addAttribute( sat::SolvAttr::name );
    q.addKind( ResKind::package );
    q.addString( pkgspec );

    MSG << "*** Get " << pkgspec << endl;
    for ( auto el : q )
      g.addNode( el );
  }

  return 0;
}


#if 0
#define msg cerr << "*** "

struct DummyStream
{
  DummyStream() : _debug( false ) {}

  template<class Tp>
  DummyStream & operator<<( Tp && val )
  { if ( _debug ) cerr << std::forward<Tp>(val); return *this; }

  DummyStream & operator<<( std::ostream& (*iomanip)( std::ostream& ) )
  { if ( _debug ) cerr << iomanip; return *this; }

  explicit operator bool() const { return _debug; }

  void turn( bool onoff_r )	{ _debug = onoff_r; }
  void on()			{ _debug = true; }
  void off()			{ _debug = false; }

  bool _debug;
} dmsg;

///////////////////////////////////////////////////////////////////

static std::string appname( "???" );

int usage( const std::string & msg_r = std::string(), int exit_r = 100 )
{
  if ( ! msg_r.empty() )
  {
    cerr << endl << msg_r << endl << endl;
  }
  cerr << "Usage: " << appname << "[-r/--root DIR] [-p] [-h]" << endl;
  cerr << "Evaporate the set of installed packages to a minimum, so " << endl;
  cerr << "their dependencies still span the whole system." << endl;
  cerr << "" << endl;
  cerr << "-h               Display this help." << endl;
  cerr << "-p               Consider only packages." << endl;
  cerr << "-r, --root DIR   Use system rooted at DIR." << endl;
  cerr << "-s, --solv FILE  Use alternate solv file." << endl;
  cerr << "-v               Verbose output." << endl;
  cerr << "" << endl;
  return exit_r;
}

int errexit( const std::string & msg_r = std::string(), int exit_r = 100 )
{
  if ( !msg_r.empty() )
    cerr << endl << msg_r << endl << endl;
  return exit_r;
}

inline bool argCheck( const char * arg, const char * value )
{ return strcmp( arg, value ) == 0; }

inline bool argCheck( const char * arg, const std::initializer_list<const char *> & values )
{ for ( const char * value : values ) if ( argCheck( arg, value ) ) return true; return false; }

#define POPARG ( ++argv,--argc )

template <typename T>
std::string whatis( T && ) { return __PRETTY_FUNCTION__; }

///////////////////////////////////////////////////////////////////
/// Simple graph with nodes representing solvables and edges representing
/// their requirements. Nodes with requirements will be chosen;  reachable
/// required nodes and their edges will be eliminated until no more
/// requirements exist.
namespace graph
{
  typedef unsigned				index_t;
  typedef std::set<index_t>			indices;
  constexpr index_t 				noindex = index_t(-1);

  bool optPackagesOnly = false;

  ///////////////////////////////////////////////////////////////////
  struct Node
  {
    Node()					: _id( noindex ) {}

    explicit operator bool() const		{ return( _id != noindex ); }
    explicit operator sat::Solvable() const	{ return sat::Solvable(_id); }

    index_t id() const				{ return _id; }

    bool inEmpty() const			{ return _in.empty(); }
    bool outEmpty() const			{ return _out.empty(); }
    size_t inCount() const			{ return _in.size(); }
    size_t outCount() const			{ return _out.size(); }
    const indices & inSet() const		{ return _in; }
    const indices & outSet() const		{ return _out; }

    bool isolated() const			{ return( inEmpty() && outEmpty() ); }

  private:
    friend class Graph;
    index_t _id;				///< Solvable id
    indices _in;				///< incoming edges
    indices _out;				///< outgoing edges
  };

  inline bool operator==( const Node & lhs, const Node & rhs )
  { return( &lhs == &rhs ); }	// actually the same instance!

  inline bool operator!=( const Node & lhs, const Node & rhs )
  { return !( lhs == rhs ); }

  std::ostream & operator<<( std::ostream & str, const Node & obj )
  {
    if ( obj )
      str << "Node[" << obj.inCount() << "->(" <<  obj.id() << ")->" << obj.outCount() << "]\t" << sat::Solvable(obj);
    else
      str << "Node[]";
    return str;
  }

  ///////////////////////////////////////////////////////////////////
  struct Graph
  {
    typedef std::vector<Node>			container_t;
    typedef container_t::const_iterator		const_iterator;
    typedef container_t::iterator		iterator;

    Graph( index_t numNodes_r  )
    : _nodes( numNodes_r + 2 )	// 0:noSolvable, 1:systemSolvable
    {
      unsigned id = unsigned(-1);
      for ( Node & node : *this )
      {
	if ( ++id >= 2 )	// kick out 0:noSolvable, 1:systemSolvableId
	  if ( !optPackagesOnly || sat::Solvable(id).isKind<Package>() )
	    node._id = id;	// Node/Solvable id is vector index
      }
    }

    bool empty() const				{ return _nodes.empty(); }
    size_t size() const				{ return _nodes.size(); }
    const_iterator begin() const		{ return _nodes.begin(); }
    const_iterator end() const			{ return _nodes.end(); }
    iterator begin()				{ return _nodes.begin(); }
    iterator end()				{ return _nodes.end(); }

    Node & operator[]( index_t id_r )		{ return _nodes[id_r]; }
    const Node & operator[]( index_t id_r ) const{ return _nodes[id_r]; }

    void addEdge( Node & from_r, Node & to_r )
    {
      // no self edge and no edges to noSolvable/systemSolvable/wiped nodes
      if ( from_r._id != to_r._id && from_r && to_r )
      {
	from_r._out.insert( to_r._id );
	to_r._in.insert( from_r._id );
      }
    }
    void addEdge( Node & from_r, index_t to_r )	{ addEdge( from_r, _nodes[to_r] ); }
    void addEdge( index_t from_r, Node & to_r )	{ addEdge( _nodes[from_r], to_r ); }
    void addEdge( index_t from_r, index_t to_r ){ addEdge( _nodes[from_r], _nodes[to_r] ); }

    /// recursively remove all outgoing nodes (protect against cycle)
    void evaporate( Node & node_r )
    {
      indices out;
      node_r._out.swap( out );	// prevent disconnecting nodes from disturbing the loop
      for ( index_t idx : out )
	rremove( _nodes[idx], node_r );
    }
    void evaporate( index_t id_r )		{ evaporate( _nodes[id_r] ); }

  public:
    // print Node indices on stream
    template <typename TStream>
    TStream & dumpIndicesOn( TStream & str, const indices & indices_r, const std::string & label_r = std::string() )
    {
      if ( ! label_r.empty() )
	str << label_r << ' ';
      str << '{';

      for ( index_t idx : indices_r )
      {
	const Node & node( _nodes[idx] );
	if ( node )
	  str << "\n  " << node;
	else
	  str << "\n  " << idx;
      }
      return str << "\n} [ indices:" << indices_r.size() << " ]";
    }

  private:
    /// recursively remove node (except protected_r)
    void rremove( Node & node_r, Node & protected_r )
    {
      if ( !node_r || node_r == protected_r )
	return;

      // disconnect from incoming nodes
      for ( index_t idx : node_r._in )
	_nodes[idx]._out.erase( node_r._id );
      node_r._in.clear();

      // recursively remove outgoing nodes
      indices out;
      node_r._out.swap( out );	// prevent disconnecting nodes from disturbing the loop
      for ( index_t idx : out )
	rremove( _nodes[idx], protected_r );

      node_r._id = noindex;	// wiped
    }

  private:
    container_t _nodes;		///< the nodes
  };

  std::ostream & operator<<( std::ostream & str, const Graph & obj )
  {
    str << "Graph {";
    size_t nodes = 0;
    size_t edges = 0;

    for ( const Node & node : obj )
    {
      if ( node )
      {
	str << "\n  " << node;
	++nodes; edges += node.outCount();
      }
      else if ( !node.isolated() )
      {
	throw("invalid and not isolated");
      }
    }
    return str << "\n} [ nodes:" << nodes << ", edges:" << edges << " ]";
  }

  ///////////////////////////////////////////////////////////////////
  /// Tarjan's strongly connected components algorithm
  namespace tarjan
  {
    struct Vertex
    {
      static constexpr unsigned novisit = unsigned(-1);

      Vertex()
      : _id( noindex )
      , _index( novisit )
      , _lowlink( novisit )
      , _onStack( false )
      {}
      Vertex( Vertex && )		= delete;	// no move - no copy
      Vertex & operator=( Vertex && )	= delete;	// no move - no copy

      explicit operator bool() const	{ return( _id != noindex ); }

      bool unvisited() const		{ return( _index == novisit ); }
      bool isRoot() const		{ return( _index == _lowlink ); }

      void visit( unsigned vidx_r )     { _index = _lowlink = vidx_r; }

      void rememberMinLowlink( unsigned seen_r )
      { if ( seen_r < _lowlink ) _lowlink = seen_r; }

      index_t		_id;		///< corresponding graph node index (or noindex)
      unsigned		_index;		///< vertex index; -1 = unprocessed
      unsigned		_lowlink;	///< smallest index of any vertex known to be reachable from this
      bool		_onStack;
    };

    std::ostream & operator<<( std::ostream & str, const Vertex & obj )
    { return str << '[' << obj._index << '.' << obj._lowlink << ']'; }

    struct Algorithm
    {
      Algorithm( const Graph & g_r, bool allsccs_r = true )
      : _g( g_r )
      , _v( g_r.size() )
      , _vidx( 0 )
      , _allsccs( allsccs_r )
      {
	// sync nodes and vertices
	for ( const Node & node : _g )
	{ if ( node ) _v[node.id()]._id = node.id(); }

	// find sccs
	for ( Vertex & v : _v )
	  if ( v && v.unvisited() )
	    scc( v );
      }

      std::list<indices> && result()	///< allow stealing the result
      { return std::move(_sccs); }

    private:
      void scc( Vertex & v )
      {
	// remember visit
	v.visit( ++_vidx );
	v._onStack = true;
	_s.push( v._id );

	// process successors; collect least Vertex::_index we can reach from here
	for ( index_t wid : _g[v._id].outSet() )
	{
	  Vertex & w( _v[wid] );
	  if ( w.unvisited() )
	  {
	    scc( w );	// depth first
	    v.rememberMinLowlink( w._lowlink );
	  }
	  else if ( w._onStack )
	  {
	    // w is on stack ==> in current scc
	    v.rememberMinLowlink( w._index );
	  }
	}

	// back at the root: collect scc from stack
	if ( v.isRoot() )
	{
	  indices scc;
	  while ( true )
	  {
	    Vertex & w( _v[_s.top()] );
	    w._onStack = false;
	    _s.pop();

	    scc.insert( w._id );

	    if ( w._id == v._id )
	      break;	// stop at v
	  }
	  if ( _allsccs || scc.size() != 1 )
	    _sccs.push_back( std::move(scc) );
	}
      }

    private:
      const Graph & _g;			///< the graph
      std::vector<Vertex> _v;		///< vertices corresponding to graph nodes
      std::stack<index_t> _s;		///< vertex stack
      unsigned _vidx;			///< current vertex visit index

      bool _allsccs;			///< whether to include single-node sccs in result
      std::list<indices> _sccs;		///< graphs list of strongly connected components
    };
  } // namespace tarjan
  ///////////////////////////////////////////////////////////////////

  /// Graphs list of strongly connected components
  std::list<indices> scc( const Graph & g, bool allsccs_r = false )
  { return tarjan::Algorithm( g, allsccs_r ).result(); }

  /// Graph as dot file
  std::ostream & writeGraphviz( std::ostream & str, const Graph & g_r, const std::string & name_r = "G" )
  {
    str << "digraph " << name_r << " {" << endl;

    for ( const Node & v : g_r )
      if ( v )
	str << v.id() << " [label=\"" << sat::Solvable(v.id()).asString() << "\"];" << endl;

    for ( const Node & v : g_r )
      if ( v && !v.outEmpty() )
	dumpRange( str << v.id() << " -> ",
		   v.outSet().begin(), v.outSet().end(),
		   "{", "", ";", "", "}" ) << endl;
    str << "}" << endl;
    return str;
  }

} // namespace graph
///////////////////////////////////////////////////////////////////

int main( int argc, const char * argv[] )
{
  appname = Pathname::basename( argv[0] );
  Pathname sysRoot("/");
  Pathname solvFile;

  while ( POPARG )
  {
    if ( argCheck( *argv, { "-h", "--help" } ) )
    {
      exit( usage( "", 0 ) );
    }
    else if ( argCheck( *argv, { "-p" } ) )
    {
      graph::optPackagesOnly = true;
    }
    if ( argCheck( *argv, { "-r", "--root" } ) )
    {
      if ( ! POPARG )
	return errexit("--root requires an argument.");
      if ( ! PathInfo( *argv ).isDir() )
	return errexit("--root requires a directory.");
      sysRoot = *argv;
    }
    else if ( argCheck( *argv, { "-s", "--solv" } ) )
    {
      if ( ! POPARG )
	return errexit("--solv requires an argument.");
      // file test delayed until --root is known
      solvFile = *argv;
    }
    else if ( argCheck( *argv, { "-v" } ) )
    {
      dmsg.on();
    }
  }

  // Go....
  sat::Pool satpool( sat::Pool::instance() );

  if ( solvFile.empty() )
  {
    // Load system at sysRoot
    msg << "Loading system at " << sysRoot << "..." << endl;
    getZYpp()->initializeTarget( sysRoot );
    getZYpp()->target()->load();
  }
  else
  {
    // Load plain solv file
    if ( sysRoot != "/" )
      solvFile = Pathname::assertprefix( sysRoot, solvFile );
    if ( ! PathInfo( solvFile ).isFile() )
      return errexit("--solv requires a file.");
    RepoInfo nrepo;
    nrepo.setAlias( sat::Pool::systemRepoAlias() );
    satpool.addRepoSolv( solvFile, nrepo );
  }

  // Build requirements graph
  msg << "Building dependency graph..." << endl;
  msg << "  Solvables: " << satpool.solvablesSize() << endl;
  unsigned dcount = 0;
  graph::Graph g( satpool.solvablesSize() );
  for ( const sat::Solvable & solv : satpool.solvables() )
  {
    graph::Node & node( g[solv.id()] );
    for ( Dep deptype : { Dep::REQUIRES , Dep::RECOMMENDS } )
    {
      for ( const Capability & dep : solv.dep( deptype ) )
      {
	sat::WhatProvides wp( dep );
	switch ( wp.size() )
	{
	  case 0:
	    if ( deptype == Dep::REQUIRES )
	      msg << "Unresolved 'requires:" << dep << "' in " << solv << endl;
	    break;
	  case 1:
	    g.addEdge( node, wp.begin()->id() );
	    ++dcount;
	    break;
	  default:
	    msg << wp.size() << " providers of '" << deptype << ":" << dep << "' in " << solv << endl;
	    break;
	}
      }
    }
  }
  graph::Graph o( g );	// save copy of the original graph
  msg << "  Dependencies: " << dcount << endl;
  dmsg << g << endl;

#if 1
  // Check for cycles
  msg << "Checking for dependency  cycles..." << endl;
  std::list<graph::indices> sccs( graph::scc( g ) );
  msg << "  Cycles: " << sccs.size() << endl;
  if ( sccs.size() )
  {
    if ( !dmsg )
      msg << "Conflating cycle of size";
    for ( const graph::indices & scc : sccs )
    {
      if ( dmsg )
	g.dumpIndicesOn( msg << "Conflating cycle of size " << scc.size(), scc, " SCC" ) << endl;
      else
	cerr << " " << scc.size() << flush;
    }
    if ( !dmsg )
      cerr << endl;
  }
#endif

  // Evaporate requirements graph
  unsigned cycle = 0;
#if 0
  for ( graph::Node & node : g )
  {
    if ( ! node.outCount() )
      continue;

    dmsg << ++cycle << " select " << node << endl;
    g.evaporate( node );
  }
#else
  while ( true )
  {
    graph::index_t id = graph::noindex;
    for ( graph::Node & node : g )
    {
      if ( node.outCount() )
      {
	if ( id == graph::noindex
	  || node.inCount() < g[id].inCount()
	  || ( node.inCount() == g[id].inCount() && node.outCount() > g[id].outCount() ) )
	  id = node.id();
      }
    }
    if ( id == graph::noindex )
      break;	// no more nodes to process

    dmsg << ++cycle << " select " << g[id] << endl;
    g.evaporate( g[id] );
  }
#endif

  {
    std::ofstream str( "G.dot" );
    graph::writeGraphviz( str, o );
  }
  {
    std::ofstream str( "E.dot" );
    graph::writeGraphviz( str, g );
  }
  msg << "Evaporated " << g << endl;
  graph::writeGraphviz( cout, g, "Evaporated" );
  return 0;
}
#endif
