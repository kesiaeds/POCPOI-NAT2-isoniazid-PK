"""Harmonize Campo Grande (CG) + Manaus (MA) into one analysis dataset.

Outputs (in /Users/kesiaeds/repos/POCPOI/analysis/data/):
  - all_participants.csv     : demographics+phenotype for ALL enrolled IDs (CG+MA)
  - pk_subset.csv            : participant-level info restricted to PK-IDs only
  - pk_nonmem_combined.csv   : NONMEM-format PK records (CG + Manaus), µg/mL

Conventions (inherited from CG INH_PGx_data.csv):
  ID prefix : CG_ or MA_ (string) to avoid numeric collisions across sites
  ACE_PROFILE : 1=Slow, 2=Intermediate, 3=Rapid
  SEX : 1=Male, 0=Female
  RACE : 1=Mixed (Parda), 2=White (Branco), 3=Black (Preta),
         4=Yellow (Amarela), 5=Indigenous (Indigena)
  TIME accumulates: Day 7 dose at t=0; Day 14 dose at t=168
  OCC : 1=Day 7 (genotype-guided), 2=Day 14 (standard 900 mg flat)
  DV in µg/mL; AMT in mg
"""
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path('/Users/kesiaeds/Library/CloudStorage/GoogleDrive-kesiaeds@stanford.edu/My Drive/Pharmacogenomics/POCPOI/Manuscript')
OUT  = Path('/Users/kesiaeds/repos/POCPOI/analysis/data')
OUT.mkdir(parents=True, exist_ok=True)

# Numeric LLOQ for isoniazid (Manaus assay reports ng/mL, LLOQ 2 ng/mL = 0.002 µg/mL).
# CG-side LLOQ per protocol/abstract is 0.01 µg/mL. We use 0.01 µg/mL globally for BLQ
# handling so the Beal M3 method / half-LLOQ substitution stays consistent across sites.
INH_LLOQ_UGML = 0.01

# Time mapping for Manaus "S2_1Hora" -> TAD in hours
PERIOD_TAD = {
    'S2_1Hora': 1,  'S2_2Hora': 2,  'S2_8Hora': 8,  'S2_24horas': 24,
    'S3_1Hora': 1,  'S3_2Hora': 2,  'S3_8Hora': 8,  'S3_24horas': 24,
}

RACE_PT_EN = {
    'Parda': 1, 'Branco': 2, 'Preta': 3, 'Indígena': 5,
    'Branca': 2, 'Amarela': 4,
}
SEX_PT = {'Masculino': 1, 'Feminino': 0}
PHENO_PT = {
    'Acetilador lento': 1,
    'Acetilador intermediário': 2,
    'Acetilador rápido': 3,
}
PHENO_EN_DEMO = {
    'Slow acetylator': 1,
    'Intermediate acetylator': 2,
    'Rapid acetylator': 3,
}
PHENO_MA_DOSE = {
    'slow': 1, 'intermediate': 2, 'intermediate*': 2, 'rapid': 3,
}

# Dose per phenotype on Day 7 / Day 14
# Day 7 (guided): slow 300, intermediate 900, rapid 1500
# Day 14: 900 flat
GUIDED_DOSE = {1: 300, 2: 900, 3: 1500}
FLAT_DOSE   = 900

# MA_3 was qPCR-classified slow, then nanopore-reclassified intermediate;
# the patient received the 300 mg slow-protocol dose on Day 7. Excluded from
# the PK analysis as a protocol deviation (genotype reclassification post-dose
# means the administered dose does not match the assigned phenotype).
EXCLUDED_PK_IDS = {'MA_3'}


# -----------------------------------------------------------------------------
# 1. Campo Grande
# -----------------------------------------------------------------------------
cg_demo = pd.read_csv(ROOT / 'CG' / 'Demographics.csv')
cg_pk   = pd.read_csv(ROOT / 'CG' / 'INH_PGx_data.csv')

cg_all = pd.DataFrame({
    'ID':           'CG_' + cg_demo['Record ID'].astype(str),
    'SITE':         'CG',
    'AGE':          cg_demo['Age'],
    'SEX':          (cg_demo['Gender'] == 'Male').astype(int),
    'RACE':         cg_demo['Race'].map({'Mixed': 1, 'White': 2, 'Black': 3, 'Yellow': 4}),
    'WT':           cg_demo['Weight'],
    'ACE_PROFILE':  cg_demo['Acetylation profile'].map(PHENO_EN_DEMO),
    'STD_DOSE':     cg_demo['Standard INH'],
    'MOD_DOSE':     cg_demo['Modified INH'],
    'HAS_PK':       cg_demo['Record ID'].isin(cg_pk['ID'].unique()),
})

# CG NONMEM PK -> prefix ID
cg_pk_nm = cg_pk.copy()
cg_pk_nm['ID']   = 'CG_' + cg_pk_nm['ID'].astype(str)
cg_pk_nm['SITE'] = 'CG'
# BLQ -> Beal M3 (matches what is done on the Manaus side above).
# DV below LLOQ becomes DV=LLOQ with CENS=1; nlmixr2 picks this up
# automatically as left-censored.
cg_pk_nm['CENS'] = 0
_blq = (cg_pk_nm['EVID'] == 0) & (cg_pk_nm['DV'] < INH_LLOQ_UGML)
cg_pk_nm.loc[_blq, 'DV']   = INH_LLOQ_UGML
cg_pk_nm.loc[_blq, 'CENS'] = 1

# -----------------------------------------------------------------------------
# 2. Manaus participants (metadata) + dose data (weight, dose plan)
# -----------------------------------------------------------------------------
ma_md   = pd.read_csv(ROOT / 'Manaus' / 'metadata_manaus.csv')
ma_dose = pd.read_excel(ROOT / 'Manaus' / 'POCPOI_dose_data_manaus.xlsx')

# Weight bug: any value > 1000 was recorded in grams -> convert to kg
ma_dose['Weight_kg'] = np.where(
    ma_dose['Weight (kg)'] > 1000,
    ma_dose['Weight (kg)'] / 1000.0,
    ma_dose['Weight (kg)']
)
ma_dose['ACE_DOSE'] = ma_dose['Acetylatorphenotype'].str.strip().str.lower().map(PHENO_MA_DOSE)

# join metadata + dose
ma_md = ma_md.rename(columns={
    'Record ID': 'rid',
    'Idade':     'AGE',
    'Gênero':    'gender_pt',
    'Raça/cor':  'race_pt',
    'Perfil de acetilação': 'pheno_pt',
})
ma_md['ACE_META'] = ma_md['pheno_pt'].map(PHENO_PT)
ma_md['SEX']      = ma_md['gender_pt'].map(SEX_PT)
ma_md['RACE']     = ma_md['race_pt'].map(RACE_PT_EN)

ma_full = ma_md.merge(
    ma_dose[['ID', 'Weight_kg', 'ACE_DOSE']].rename(columns={'ID': 'rid'}),
    on='rid', how='left'
)
# Phenotype: prefer the corrected (dose-data) value (ID 3 reclassified slow->intermediate)
ma_full['ACE_PROFILE'] = ma_full['ACE_DOSE'].fillna(ma_full['ACE_META']).astype('Int64')

# PK IDs
pk_xl = pd.ExcelFile(ROOT / 'Manaus' / 'Manaus_INH_RPT_ac-INH_Clinical_Samples_preliminary.xlsx')
s2 = pd.read_excel(pk_xl, 'S2_Summary')
s3 = pd.read_excel(pk_xl, 'S3_Summary')

def _clean_pk(df):
    d = df[['record_id', 'id', 'period', 'Isoniazid (ng/mL)']].copy()
    d['id'] = pd.to_numeric(d['id'], errors='coerce')
    d = d[d['id'].notna()].copy()
    d['id'] = d['id'].astype(int)
    d['DV_ngml'] = pd.to_numeric(d['Isoniazid (ng/mL)'], errors='coerce')
    d = d[d['period'].isin(PERIOD_TAD)]
    return d[['id', 'period', 'DV_ngml']]

ma_s2 = _clean_pk(s2)
ma_s3 = _clean_pk(s3)
ma_pk_ids = sorted(set(ma_s2['id']) | set(ma_s3['id']))

ma_all = pd.DataFrame({
    'ID':           'MA_' + ma_full['rid'].astype(str),
    'SITE':         'MA',
    'AGE':          ma_full['AGE'],
    'SEX':          ma_full['SEX'].astype('Int64'),
    'RACE':         ma_full['RACE'].astype('Int64'),
    'WT':           ma_full['Weight_kg'],
    'ACE_PROFILE':  ma_full['ACE_PROFILE'],
    'STD_DOSE':     FLAT_DOSE,
    'MOD_DOSE':     ma_full['ACE_PROFILE'].map(GUIDED_DOSE).astype('Int64'),
    'HAS_PK':       ma_full['rid'].isin(ma_pk_ids),
})
# Drop excluded participants from the PK subset (kept in enrolled cohort).
ma_all.loc[ma_all['ID'].isin(EXCLUDED_PK_IDS), 'HAS_PK'] = False

# -----------------------------------------------------------------------------
# 3. Combined participant table (ALL IDs)
# -----------------------------------------------------------------------------
all_participants = pd.concat([cg_all, ma_all], ignore_index=True)
all_participants.to_csv(OUT / 'all_participants.csv', index=False)

# PK-only subset
pk_subset = all_participants[all_participants['HAS_PK']].copy()
pk_subset.to_csv(OUT / 'pk_subset.csv', index=False)

# -----------------------------------------------------------------------------
# 4. Build Manaus NONMEM-format PK records (TIME accumulating, OCC 1/2)
# -----------------------------------------------------------------------------
records = []
ma_dose_lookup = ma_all.set_index('ID')[['WT', 'AGE', 'RACE', 'SEX', 'ACE_PROFILE', 'MOD_DOSE']]

def _emit(rid, day, occ, dose, obs_rows, demo):
    """Emit one occasion: a dose row at t_base, then observation rows."""
    t_base = 0 if day == 7 else 168
    # dose record
    records.append({
        'ID': rid, 'TIME': t_base, 'DV': 0.0, 'AMT': int(dose), 'MDV': 1, 'EVID': 1,
        'CENS': 0,
        'TAD': 0, 'DAY': day, 'DOSE': int(dose / demo['WT']) if demo['WT'] else np.nan,
        'ACE_PROFILE': int(demo['ACE_PROFILE']), 'WT': float(demo['WT']),
        'AGE': int(demo['AGE']),
        'RACE': int(demo['RACE']) if pd.notna(demo['RACE']) else np.nan,
        'SEX':  int(demo['SEX'])  if pd.notna(demo['SEX'])  else np.nan,
        'OCC':  occ, 'SITE': 'MA',
    })
    for _, r in obs_rows.iterrows():
        tad = PERIOD_TAD[r['period']]
        dv_ngml = r['DV_ngml']
        if pd.isna(dv_ngml):
            continue  # drop missing
        dv_ugml = dv_ngml / 1000.0
        mdv = 0
        # BLQ -> Beal M3: set DV to LLOQ, mark CENS=1; nlmixr2 picks this up
        # automatically and uses the left-censored likelihood instead of the
        # half-LLOQ substitution.
        if dv_ugml < INH_LLOQ_UGML:
            dv_ugml = INH_LLOQ_UGML
            cens = 1
        else:
            cens = 0
        records.append({
            'ID': rid, 'TIME': t_base + tad, 'DV': round(dv_ugml, 4), 'AMT': 0,
            'MDV': mdv, 'EVID': 0, 'CENS': cens, 'TAD': tad, 'DAY': day,
            'DOSE': int(dose / demo['WT']) if demo['WT'] else np.nan,
            'ACE_PROFILE': int(demo['ACE_PROFILE']), 'WT': float(demo['WT']),
            'AGE': int(demo['AGE']),
            'RACE': int(demo['RACE']) if pd.notna(demo['RACE']) else np.nan,
            'SEX':  int(demo['SEX'])  if pd.notna(demo['SEX'])  else np.nan,
            'OCC':  occ, 'SITE': 'MA',
        })

for rid_num in ma_pk_ids:
    rid = f'MA_{rid_num}'
    if rid not in ma_dose_lookup.index:
        continue  # no dose/demo info, skip
    demo = ma_dose_lookup.loc[rid]
    if pd.isna(demo['ACE_PROFILE']) or pd.isna(demo['WT']):
        continue
    # Day 7 = genotype-guided
    guided = GUIDED_DOSE[int(demo['ACE_PROFILE'])]
    s2_rows = ma_s2[ma_s2['id'] == rid_num].sort_values('period')
    if len(s2_rows):
        _emit(rid, 7, 1, guided, s2_rows, demo)
    # Day 14 = 900 mg flat
    s3_rows = ma_s3[ma_s3['id'] == rid_num].sort_values('period')
    if len(s3_rows):
        _emit(rid, 14, 2, FLAT_DOSE, s3_rows, demo)

ma_pk_nm = pd.DataFrame(records, columns=cg_pk_nm.columns)

# -----------------------------------------------------------------------------
# 5. Combined NONMEM dataset
# -----------------------------------------------------------------------------
combined = pd.concat([cg_pk_nm, ma_pk_nm], ignore_index=True)
combined = combined[~combined['ID'].isin(EXCLUDED_PK_IDS)].copy()
combined = combined.sort_values(['ID', 'TIME', 'EVID'], kind='stable').reset_index(drop=True)
combined.to_csv(OUT / 'pk_nonmem_combined.csv', index=False)

# -----------------------------------------------------------------------------
# 6. Summary
# -----------------------------------------------------------------------------
def _phen_counts(df):
    s = df['ACE_PROFILE'].map({1: 'slow', 2: 'intermediate', 3: 'rapid'}).value_counts()
    return s.reindex(['slow', 'intermediate', 'rapid'], fill_value=0).to_dict()

print('=== ALL participants ===')
print('Total:', len(all_participants), 'CG:', (all_participants['SITE']=='CG').sum(),
      'MA:', (all_participants['SITE']=='MA').sum())
print('Phenotype (all):', _phen_counts(all_participants))
print('Phenotype CG:',    _phen_counts(all_participants[all_participants['SITE']=='CG']))
print('Phenotype MA:',    _phen_counts(all_participants[all_participants['SITE']=='MA']))
print()
print('=== PK subset ===')
print('Total:', len(pk_subset), 'CG:', (pk_subset['SITE']=='CG').sum(),
      'MA:', (pk_subset['SITE']=='MA').sum())
print('Phenotype:', _phen_counts(pk_subset))
print()
print('=== Combined NONMEM dataset ===')
print('Rows:', len(combined), 'unique IDs:', combined['ID'].nunique())
print('Observations (EVID=0):', (combined['EVID']==0).sum())
print('Dose events (EVID=1):',  (combined['EVID']==1).sum())
print()
print('Files written to:', OUT)
