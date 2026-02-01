import pandas as pd

# Archivos
file_total = 'tproduct_total.out'
file_tail  = 'tproduct_tail.out'

# Regiones en el orden del tally
regions = [110, 111, 112, 200, 230, 300]

def read_sum_over(file_name, regions):
    values = []
    with open(file_name, "r") as f:
        for line in f:
            if line.strip().startswith("#   sum over"):
                values.append(float(line.split()[-1]))

    if len(values) != len(regions):
        raise ValueError(
            f"{file_name}: encontrados {len(values)} sum over, "
            f"pero se esperaban {len(regions)}"
        )

    return pd.DataFrame({
        "reg": regions,
        "N_per_source": values
    })

# Leer ambos archivos
df_total = read_sum_over(file_total, regions)
df_tail  = read_sum_over(file_tail, regions)

# Combinar y calcular fracción
df = df_total.merge(df_tail, on="reg", suffixes=("_total", "_tail"))
df["tail_fraction"] = df["N_per_source_tail"] / df["N_per_source_total"]

print(df)

