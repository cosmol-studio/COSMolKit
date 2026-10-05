import cosmolkit as ck
import pytest

CASES = [
    (ck.AlignmentParameters,dict(probe_conformer_id=-1,reference_conformer_id=-1,atom_map=None,weights=None,reflect=False,max_iterations=50)),
    (ck.BestAlignmentParameters,dict(probe_conformer_id=-1,reference_conformer_id=-1,atom_maps=[],weights=None,reflect=False,max_iterations=50,max_matches=1000000,symmetrize_conjugated_terminal_groups=True,ignore_hydrogens=True,num_threads=1)),
    (ck.CoordinateRmsdParameters,dict(probe_conformer_id=-1,reference_conformer_id=-1,atom_maps=[],weights=None,max_matches=1000000,symmetrize_conjugated_terminal_groups=True)),
    (ck.AllConformerRmsdParameters,dict(atom_maps=[],weights=None,max_matches=1000000,symmetrize_conjugated_terminal_groups=True,ignore_hydrogens=True,num_threads=1)),
    (ck.ConformerAlignmentParameters,dict(atom_indices=None,conformer_ids=None,weights=None,reflect=False,max_iterations=50)),
]
@pytest.mark.parametrize("typ,defaults",CASES)
def test_all_original_parameter_fields_defaults_and_mutable_state(typ,defaults):
    p=typ(); assert type(p) is typ
    for k,v in defaults.items():assert getattr(p,k)==v
    edits=dict(defaults)
    for k in edits:
        if k in ["probe_conformer_id","reference_conformer_id"]:edits[k]=17
        elif k=="num_threads":edits[k]=-1
        elif k=="max_matches":edits[k]=0
        elif k=="max_iterations":edits[k]=0
        elif k=="weights":edits[k]=[1.,2.]
        elif k in ["atom_indices","conformer_ids"]:edits[k]=[0,17]
        elif k=="atom_map":edits[k]=[ck.AlignmentAtomMap(0,1),ck.AlignmentAtomMap(1,0)]
        elif k=="atom_maps":edits[k]=[[ck.AlignmentAtomMap(0,1),ck.AlignmentAtomMap(1,0)]]
        elif isinstance(edits[k],bool):edits[k]=not edits[k]
    q=typ(**edits)
    for k,v in edits.items():
        setattr(p,k,v)
        def scalar(x):
            if isinstance(x,ck.AlignmentAtomMap):return (x.probe_atom,x.reference_atom)
            if isinstance(x,list):return [scalar(a) for a in x]
            return x
        assert scalar(getattr(p,k))==scalar(getattr(q,k))==scalar(v)
    for k in ["weights","atom_indices","conformer_ids"]:
        if k in edits:
            detached=getattr(p,k);detached.append(999);assert getattr(p,k)==edits[k]

def molecule():
    b=ck.Molecule.from_smiles("CC").to_builder();b.add_3d_conformer([[0.,0.,0.],[1.,0.,0.]]);b.add_3d_conformer([[2.,1.,0.],[3.,1.,0.]]);return b.build()

def test_original_atom_map_constructor_repr_and_fields():
    p=ck.AlignmentAtomMap(2,7);assert repr(p)=="AlignmentAtomMap(probe_atom=2, reference_atom=7)";assert (p.probe_atom,p.reference_atom)==(2,7)
    # Original pyclass uses get_all/set_all, including atom-map fields.
    p.probe_atom=3;p.reference_atom=8;assert (p.probe_atom,p.reference_atom)==(3,8)
    with pytest.raises(OverflowError):ck.AlignmentAtomMap(-1,0)

def test_all_original_result_protocols_and_detached_reads():
    m=molecule();r=m.alignment_transform_to(m);assert isinstance(r,ck.AlignmentResult);assert repr(r)=="AlignmentResult(rmsd=0, mapped_atoms=2)";assert r.rmsd()==0.
    t=r.transform();assert isinstance(t,ck.AlignmentTransform);assert repr(t)=="AlignmentTransform(matrix=4x4)";a=t.matrix();assert len(a)==4 and all(len(x)==4 for x in a);a[0][0]=99;assert t.matrix()[0][0]==1.
    a=r.atom_map();assert len(a)==2;a.clear();assert len(r.atom_map())==2
    rows=m.all_conformer_best_rmsds();assert len(rows)==1;v=rows[0];assert isinstance(v,ck.ConformerRmsd);assert repr(v)=="ConformerRmsd(probe_conformer_id=1, reference_conformer_id=0, rmsd=0)";assert (v.probe_conformer_id(),v.reference_conformer_id(),v.rmsd())==(1,0,0.)
    aligned,report=m.with_aligned_conformers();assert isinstance(report,ck.ConformerAlignmentReport);assert repr(report)=="ConformerAlignmentReport(rmsds=1)";a=report.rmsds();assert a==[0.];a[0]=99;assert report.rmsds()==[0.]
    for typ in [ck.AlignmentResult,ck.AlignmentTransform,ck.ConformerRmsd,ck.ConformerAlignmentReport]:
        with pytest.raises(TypeError):typ()

def test_configured_public_forms_use_same_report_and_typed_errors():
    m=molecule();ref=molecule();p=ck.AlignmentParameters();bp=ck.BestAlignmentParameters();cp=ck.CoordinateRmsdParameters();ap=ck.AllConformerRmsdParameters();fp=ck.ConformerAlignmentParameters()
    assert m.alignment_transform_to_with_params(ref,p).rmsd()==m.alignment_transform_to(ref).rmsd()
    assert m.best_alignment_to_with_params(ref,bp).rmsd()==m.best_alignment_to(ref).rmsd()
    assert m.best_rmsd_to_with_params(ref,bp)==m.best_rmsd_to(ref)
    assert m.coordinate_rmsd_to_with_params(ref,cp)==m.coordinate_rmsd_to(ref)
    assert len(m.all_conformer_best_rmsds_with_params(ap))==len(m.all_conformer_best_rmsds())
    assert m.with_alignment_to_with_params(ref,p)[1].rmsd()==0.
    assert m.with_aligned_conformers_with_params(fp)[1].rmsds()==[0.]
    assert m.align_to_with_params_(ref,p).rmsd()==0.
    assert m.align_conformers_with_params_(fp).rmsds()==[0.]
    bad=ck.AlignmentParameters(atom_map=[ck.AlignmentAtomMap(7,0)])
    before=m.to_binary()
    with pytest.raises(ck.AlignmentError) as e:m.alignment_transform_to(ref,bad)
    assert e.value.domain=="alignment" and e.value.kind=="ProbeAtomOutOfRange" and e.value.index==7 and e.value.atom_count==2
    with pytest.raises(ck.OperationError) as e:m.align_to_(ref,bad)
    assert e.value.kind=="Alignment" and isinstance(e.value.__cause__,ck.AlignmentError)
    assert e.value.__cause__.index==7 and m.to_binary()==before
