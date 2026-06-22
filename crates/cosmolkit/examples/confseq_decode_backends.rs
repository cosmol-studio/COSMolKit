use cosmolkit::{ConfSeqDecodeOptions, ConfSeqTemplateBackend, decode_confseq_record_with_options};

// Tokenized corpus record: dude_aa2ar__L021087 candidate 0.
// Spaces are ConfSeq token boundaries and must be preserved.
const CONFSEQ: &str = "N C <112> | ( = O ) <84> C <111> | <-45> n 1 c c ( <21> N <123> | <173> C <116> | ( = O ) <120> c 2 c c c c c 2 <112> C <112> | <-66> N 2 <-172> C ( = O ) <-2> C <0> N <3> C <174> 2 = O ) c n 1";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let reference = decode_confseq_record_with_options(
        CONFSEQ,
        &ConfSeqDecodeOptions {
            optimize_with_uff: false,
            template_backend: ConfSeqTemplateBackend::DistanceGeometry,
            ..ConfSeqDecodeOptions::default()
        },
    )?;

    let fast = decode_confseq_record_with_options(
        CONFSEQ,
        &ConfSeqDecodeOptions {
            optimize_with_uff: false,
            template_backend: ConfSeqTemplateBackend::FastGeometry,
            ..ConfSeqDecodeOptions::default()
        },
    )?;

    println!("Corpus record: dude_aa2ar__L021087 candidate 0");
    println!("ConfSeq numeric tokens: {}", CONFSEQ.matches('<').count());
    println!(
        "DistanceGeometry decoded atoms={} conformers={}",
        reference.num_atoms(),
        reference.conformers_3d().len()
    );
    println!(
        "FastGeometry decoded atoms={} conformers={}",
        fast.num_atoms(),
        fast.conformers_3d().len()
    );

    Ok(())
}
