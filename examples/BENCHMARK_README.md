# Benchmark Scripts

Este directorio contiene scripts para evaluar el rendimiento de sunwhere comparado con otras bibliotecas de cálculo de posición solar.

## Scripts Disponibles

### 1. benchmark_simple.py

Script básico que solo requiere sunwhere. Evalúa el rendimiento interno de sunwhere con diferentes algoritmos y números de sitios.

**Uso:**
```bash
python benchmark_simple.py
```

**Características:**
- ✓ No requiere bibliotecas externas (solo sunwhere)
- ✓ Prueba 1 y 100 sitios
- ✓ Algoritmos PSA y NREL
- ✓ Motor numexpr
- ✓ Métricas de rendimiento y throughput

**Salida esperada:**
- Tiempos de ejecución para cada configuración
- Comparación de eficiencia de vectorización
- Throughput (cálculos por segundo)

---

### 2. benchmark_comparison.py

Script completo que compara sunwhere con pvlib y solpox.

**Requisitos:**
```bash
pip install sunwhere pvlib solpox tabulate
```

**Uso:**
```bash
python benchmark_comparison.py
```

**Características:**
- ✓ Compara 3 bibliotecas: sunwhere, pvlib, solpox
- ✓ Prueba 1 y 100 sitios
- ✓ Algoritmos PSA y NREL
- ✓ Tabla comparativa de resultados
- ✓ Análisis de speedup relativo

**Notas:**
- pvlib: Optimizada para un solo sitio, no vectoriza bien sobre múltiples sitios
- solpox: Puede no soportar múltiples sitios
- sunwhere: Optimizada específicamente para múltiples sitios

---

## Configuración de Pruebas

Ambos scripts usan la misma configuración base:

- **Periodo de tiempo:** 1 año (8760 timestamps horarios)
- **Sitios:** 1 y 100 ubicaciones
- **Algoritmos:** PSA y NREL
- **Motor:** numexpr (cuando está disponible)
- **Repeticiones:** 5 ejecuciones por configuración
- **Sin refracción:** Para comparación más justa

## Interpretación de Resultados

### Métricas Clave

1. **Tiempo medio (ms):** Tiempo promedio de ejecución
2. **Desviación estándar (ms):** Variabilidad en las mediciones
3. **Speedup:** Factor de mejora respecto a ejecución secuencial
4. **Eficiencia:** Porcentaje de vectorización efectiva
5. **Throughput:** Cálculos por segundo

### Ejemplo de Resultados Esperados

Para 100 sitios con PSA:
- **sunwhere (numexpr):** ~50-100 ms (muy eficiente)
- **pvlib:** No optimizado para múltiples sitios
- **Eficiencia de vectorización:** >80% es excelente

## Ejecutar con el Entorno Virtual

Si estás en el directorio del proyecto:

```bash
# Usando uv
uv run python benchmark_simple.py

# O directamente con el venv
.venv/bin/python benchmark_simple.py
```

## Troubleshooting

### Error: "No module named 'pvlib'"
```bash
pip install pvlib
```

### Error: "No module named 'solpox'"
```bash
pip install solpox
```

### Error: "No module named 'tabulate'"
```bash
pip install tabulate
```
O ejecuta `benchmark_simple.py` que no requiere tabulate.

## Benchmarks Adicionales

El proyecto también incluye benchmarks más completos en:
- `src/sunwhere/_cli/benchmark.py` - Benchmarks contra múltiples bibliotecas incluyendo SolTrack, sg2, SPARTA

Para ejecutar los benchmarks CLI:
```bash
sunwhere benchmark --help
```
