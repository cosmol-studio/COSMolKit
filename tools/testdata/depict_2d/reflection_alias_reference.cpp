// Fixed native reference only. RDKit pin351f8f378f8ad6bbd517980c38896e66bf907af8c.
// Calls original EmbeddedAtom::Reflect and EmbeddedFrag::flipAboutBond.
// No geometry/caller algorithm is implemented here.
#include <GraphMol/Depictor/EmbeddedFrag.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <array>
#include <cstdint>
#include <cstring>
#include <iostream>

using RDDepict::EmbeddedAtom;
using RDDepict::EmbeddedFrag;
using RDGeom::Point2D;
std::uint64_t bits(double value) {
  std::uint64_t result;
  std::memcpy(&result, &value, sizeof(result));
  return result;
}
double value(std::uint64_t bits) {
  double result;
  std::memcpy(&result, &bits, sizeof(result));
  return result;
}
void row(const char *family, unsigned id, const char *stage, unsigned key,
         const EmbeddedAtom &atom) {
  std::cout << family << '\t' << id << '\t' << stage << '\t' << key << '\t'
    << atom.aid << '\t' << bits(atom.loc.x) << '\t' << bits(atom.loc.y) << '\t'
    << bits(atom.normal.x) << '\t' << bits(atom.normal.y) << '\t' << atom.ccw
    << '\t' << bits(atom.angle) << '\t' << atom.nbr1 << '\t' << atom.nbr2
    << '\t' << atom.CisTransNbr << '\t' << atom.rotDir << '\t'
    << bits(atom.d_density) << '\t' << atom.df_fixed << '\t';
  for (unsigned i=0;i<atom.neighs.size();++i) {
    if (i) std::cout << ',';
    std::cout << atom.neighs[i];
  }
  std::cout << '\n';
}
void kernel(unsigned id, Point2D first, Point2D second, Point2D distinct,
            Point2D normal, unsigned alias, bool ccw) {
  EmbeddedAtom atom;
  atom.loc = alias == 1 ? first : alias == 2 ? second : distinct;
  atom.normal = normal;
  atom.ccw = ccw;
  row("K", id, "PRE", 0, atom);
  const Point2D &a = alias == 1 ? atom.loc : first;
  const Point2D &b = alias == 2 ? atom.loc : second;
  atom.Reflect(a, b);
  row("K", id, "POST", 0, atom);
}
int main() {
  unsigned calls=0;
  for (unsigned frame=0;frame<2;++frame) {
    const Point2D a = frame ? Point2D(10.17,-3.41) : Point2D(-.23,1.17);
    const Point2D b = frame ? Point2D(12.51,2.87) : Point2D(2.51,-.87);
    const Point2D loc = frame ? Point2D(11.13,4.29) : Point2D(3.14,1.59);
    const Point2D normal = frame ? Point2D(-.37,.91) : Point2D(.61,.43);
    for (unsigned alias=0;alias<3;++alias) for (bool ccw : {false,true})
      kernel(calls++,a,b,loc,normal,alias,ccw);
  }
  const std::uint64_t retained[2][9] = {
    {4624945975734500932ULL,4618271276495847798ULL,4606722574913418428ULL,13822726034664128445ULL,4624679614270839367ULL,4616668647364321170ULL,4624945975734500932ULL,4618271276495847798ULL,4607182418800017408ULL},
    {4624945975734500934ULL,4618271276495847797ULL,13830094611768194224ULL,4599353997809352688ULL,4624679614270839367ULL,4616668647364321170ULL,4624945975734500934ULL,4618271276495847797ULL,0ULL}
  };
  for (const auto &r : retained)
    kernel(calls++,Point2D(value(r[4]),value(r[5])),Point2D(value(r[6]),value(r[7])),
           Point2D(value(r[0]),value(r[1])),Point2D(value(r[2]),value(r[3])),2,value(r[8])!=0);
  const Point2D locs[] = {{-.23,1.17},{2.51,-.87},{3.14,1.59},{-2.73,.41},{1.13,2.29}};
  const unsigned edges[2][4][2] = {{{0,1},{1,2},{2,3},{3,4}},{{0,1},{1,2},{1,3},{3,4}}};
  unsigned caller_calls=0;
  for (unsigned tree=0;tree<2;++tree) for (bool flip : {false,true})
    for (unsigned fixed=0;fixed<3;++fixed) {
      RDKit::RWMol mol;
      for (unsigned i=0;i<5;++i) mol.addAtom(new RDKit::Atom(6),true,true);
      for (const auto &edge : edges[tree]) mol.addBond(edge[0],edge[1],RDKit::Bond::SINGLE);
      mol.updatePropertyCache(false);
      RDKit::MolOps::fastFindRings(mol);
      RDGeom::INT_POINT2D_MAP coords;
      for (unsigned i=0;i<5;++i) coords[i]=locs[i];
      EmbeddedFrag fragment(&mol,coords);
      // Legal fixture mutation: original object is NONCONST; accessor is const.
      auto &atoms=const_cast<RDDepict::INT_EATOM_MAP &>(fragment.GetEmbeddedAtoms());
      for (auto &[key,atom] : atoms) {
        atom.aid=0; atom.loc=locs[key]; atom.normal=Point2D(.61,.43);
        atom.ccw=key%2; atom.angle=-1.; atom.nbr1=-1; atom.nbr2=-1;
        atom.CisTransNbr=-1; atom.rotDir=0; atom.neighs.clear();
        atom.d_density=-1.; atom.df_fixed=fixed && key==fixed;
        row("C",caller_calls,"PRE",key,atom);
      }
      fragment.flipAboutBond(1,flip);
      for (const auto &[key,atom] : fragment.GetEmbeddedAtoms())
        row("C",caller_calls,"POST",key,atom);
      ++caller_calls;
    }
  std::cout << "CENSUS\t" << calls << '\t' << caller_calls << '\n';
  return calls==14 && caller_calls==12 ? 0 : 1;
}
