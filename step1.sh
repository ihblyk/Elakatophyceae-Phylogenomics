python - <<'PY'
import pandas as pd
f="/mnt/sdb/transcriptome/elakatophyceae_pep/pp/deredundancy/0212/cross/OrthoFinder/Orthogroups/Orthogroups.tsv"
df=pd.read_csv(f,sep="\t",index_col=0)

def ncopy(cell):
    if pd.isna(cell) or str(cell).strip()=="": return 0
    return len([x for x in str(cell).split(", ") if x.strip()])

def pick(min_sp,max_cp,name):
    keep=[]
    for og,row in df.iterrows():
        cp=[ncopy(row[c]) for c in df.columns]
        present=sum(x>0 for x in cp)
        if present>=min_sp and all(x<=max_cp for x in cp):
            keep.append(og)
    out=f"/mnt/sdb/transcriptome/elakatophyceae_pep/pp/deredundancy/0212/cross/OrthoFinder/Orthogroups/{name}.txt"
    open(out,"w").write("\n".join(keep))
    print(name,len(keep),out)

pick(5,1,"SOG_5of5_1copy")
pick(4,1,"SOG_4of5_1copy")
PY