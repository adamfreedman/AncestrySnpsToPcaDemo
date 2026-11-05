# Trait Loci and Genotypes: Overview

This dataset reports genotypes at 21 trait loci for each dog (`ccny-genomics_selected_trait_genotypes.csv`). The meaning of each genotype is in `trait_genotype_summary.csv`, one row per genotype, with the columns `trait_gene_id`, `genotype`, `meaning` and `gene_interactions`.

## Loci by category

**Coat color and pattern**

| Locus | What it influences |
|---|---|
| `MC1R_E_locus` | Whether eumelanin pigment (black/brown) can be made in the fur. Also controls masks and similar patterns. This locus alters how most other color loci are expressed. |
| `CBD103_K_locus` | Solid dark coat (K<sup>B</sup>) or a patterned coat (k<sup>y</sup>) |
| `ASIP_A_locus` | Pattern type: fawn sable, agouti (wolf sable), black/brown and tan, or recessive black/brown |
| `RALY_Saddle_trait_gene` | Saddle tan, a variant of the tan-point pattern |
| `TYRP1_B_locus` | Black or brown pigment |
| `MLPH_D_locus` | Dilution, which lightens black to gray and brown to a lighter brown |
| `Intensity_red_pigment` | Shade of red/yellow pigment, from red through tan to cream |
| `PMEL_Merle` | Merle and double merle patterning |
| `MITF_white_spotting` | Amount of white spotting |
| `USH2A_Roan` | Roan (flecks of color within white areas) |

**Coat type:** `FGF5_longcoat` (coat length), `KRT71_CurlyCoat` (curl), `RSPO2_moustache` (furnishings: beard, mustache, eyebrows, wirehair)

**Body size:** `GHR_size1`, `GHR_size2`, `IGF1_size`, `IGFR1_toy`, `STC2_size`. Each genotype is labeled *Smaller*, *Intermediate* or *Larger*. Body size is polygenic and the prediction you see in your dog's profile comes from a model using hundreds of markers. These loci here are some with the largest effect sizes.

**Other:** `ALX4_Blue_Eyes_Linkage` (blue eyes), `LMBR1_Claw` (hind dewclaws), `EPAS1_altitude` (tolerance to high altitude)

## Reading the genotypes

- Two alleles are listed per genotype, e.g. `AyAt` or `Ssp`.
- `X or Y` (e.g. `Bb or bb`, `Eme or Ee`) marks a call that could not be fully resolved.
- In `PMEL_Merle`, `M*` stands for any Merle allele. The Merle phenotype results from SINE insertions into the PMEL gene. Multiple length polymorphisms have been observed.
- A blank or `No Call` means no genotype is available for that dog at that locus.
- `MC1R_E_locus` genotypes have no standalone meaning. Their effect shows up through the other coat color loci.

## Gene interactions

Coat color comes from several loci acting in an epistatic hierarchy. For example, an `ee` genotype at the E-locus stops dark pigment from forming in the fur. This hides the effects of the K-, A- and B-loci and the merle locus in the fur. Likewise, a K<sup>B</sup> allele hides the A-locus pattern. The `gene_interactions` column says which other loci may change how each genotype is expressed, but it does not give the full logic.

For a full explanation of how the coat color loci interact, see Embark's [Coat Color Genetics 101](https://embarkvet.com/blogs/resources/science-corner-coat-color-genetics-101).
