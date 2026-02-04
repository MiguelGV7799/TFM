# PHITSFusion © 2026 by Miguel García, Manuel Cotelo
# Licensed under CC BY-NC-SA 4.0
# To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/


import pandas as pd

# Archivos
file_atomic = 'tproduct_photon_atomic.out'
file_nuclear  = 'tproduct_photon_nuclear.out'

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
df_atomic  = read_sum_over(file_atomic, regions).rename(
    columns={"N_per_source": "photons_atomic"}
)
df_nuclear = read_sum_over(file_nuclear, regions).rename(
    columns={"N_per_source": "photons_nuclear"}
)
# Combinar por región
df = df_atomic.merge(df_nuclear, on="reg")

# Suma total por región
df["photons_total"] = df["photons_atomic"] + df["photons_nuclear"]

# Fila de totales globales
total_row = pd.DataFrame([{
    "reg": "TOTAL",
    "photons_atomic": df["photons_atomic"].sum(),
    "photons_nuclear": df["photons_nuclear"].sum(),
    "photons_total": df["photons_total"].sum(),
}])

# Tabla final
df_out = pd.concat([df, total_row], ignore_index=True)

# Mostrar con formato científico
pd.set_option("display.float_format", "{:.6e}".format)

print(df_out)
