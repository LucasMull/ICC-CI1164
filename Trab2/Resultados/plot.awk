#!/usr/bin/awk -f

BEGIN {
  FS = ","
  # size: 10, 32, 50, 64, 100, 128, 200, 256, 300, 400, 512, 1000
  # type: FLOPS_DP, L2CACHE, L3
  # table: AjusteDeCurvas; Triangulariza; TriangularizaOtimiz
  # op $1: valor especifico; op 1: runtime counter
  graph[size][type][table][op] = ""
}

BEGINFILE {
  type = ""
 sdhasqw
 size = ""
  table = ""

  match(FILENAME, /_[0-9]*.txt$/)
  switch (type = substr(FILENAME, 0, RSTART-1)) {
  case "FLOPS_DP":
  case "L2CACHE":
  case "L3":
      size = substr(FILENAME, RSTART+1, RLENGTH-5)
      break
  default:
      print FILENAME" será ignorado" > "/dev/stderr"
      nextfile
  }
}

FNR > 6 { # pula primeiras 6 linhas de cabeçalho
  switch ($1) {
  case "TABLE":
      table = $2
      next
  case "L3 bandwidth [MBytes/s]":
  case "L2 miss ratio":
  case "DP MFLOP/s":
  case "AVX DP MFLOP/s":
      graph[size][type][table][$1] = $2
      next
  case "RDTSC Runtime [s]":
      if (table == "Region Interpolacao") 
        next

      graph[size][type][table][$1] += $2
      graph[size][type][table][1]++

      next
  default:
      next
  }
}

ENDFILE { 
  # Calcula media dos tempos
  op = "RDTSC Runtime [s]"
  if (graph[size][type][table][op]) {
    graph[size][type][table][op] = graph[size][type][table][op] / graph[size][type][table][1]
  }
}

END { 

  graph[size][type][table][op]

  # Imprime cabeçalho Run-Time AjusteDeCurvas
  print "# N AjusteDeCurvas" > "tempo_ajc.csv"
  for (isize in graph) {
    if (!isize) continue

    printf "%s ", isize >> "tempo_ajc.csv"

    runtime = 0
    for (itype in graph[isize]) {
      if (!itype) continue

      runtime += graph[isize][itype]["Region AjusteDeCurvas"]["RDTSC Runtime [s]"];
    }

    if (runtime)
      printf "%f\n", runtime / 3 >> "tempo_ajc.csv"
  }

  # Imprime cabeçalho para Banda de Memória
  print "# MemoryBandwidth[MBytes/s]" > "mem_ajc.csv"
  print "# N AjusteDeCurvas" >> "mem_ajc.csv"
  for (isize in graph) {
      if (!isize) continue

      printf "%s ", isize >> "mem_ajc.csv"

      bandwidth = graph[isize]["L3"]["Region AjusteDeCurvas"]["L3 bandwidth [MBytes/s]"]

      printf "%f\n", bandwidth >> "mem_ajc.csv"
  }

  # Imprime cabeçalho para Cache miss L2
  print "# DataCacheMissRatio" > "cmiss_ajc.csv"
  print "# N AjusteDeCurvas" >> "cmiss_ajc.csv"
  for (isize in graph) {
      if (!isize) continue

      printf "%s ", isize >> "cmiss_ajc.csv"

      cmiss = graph[isize]["L2CACHE"]["Region AjusteDeCurvas"]["L2 miss ratio"]

      printf "%f\n", cmiss >> "cmiss_ajc.csv"
  }

  # Imprime cabeçalho para Cache miss L2
  print "#   FLOPS_DP" > "flops_ajc.csv"
  print "# N AjusteDeCurvas" >> "flops_ajc.csv"
  for (isize in graph) {
      if (!isize) continue

      printf "%s ", isize >> "flops_ajc.csv"

      flops = graph[isize]["FLOPS_DP"]["Region AjusteDeCurvas"]["DP MFLOP/s"]

      printf "%f\n", flops >> "flops_ajc.csv"
  }

  # Imprime cabeçalho Run-Time Triangulariza e TriangularizaOtimiz
  print "# N TriangularizaOtimiz Triangulariza" > "tempo_triang.csv"
  for (isize in graph) {
    if (!isize) continue

    printf "%s ", isize >> "tempo_triang.csv"

    runtime_otimiz = 0
    runtime_normal = 0
    for (itype in graph[isize]) {
      if (!itype) continue

      runtime_otimiz += graph[isize][itype]["Region TriangularizaOtimiz"]["RDTSC Runtime [s]"];
      runtime_normal += graph[isize][itype]["Region Triangulariza"]["RDTSC Runtime [s]"];
    }

    if (runtime_otimiz && runtime_normal)
      printf "%f %f\n", runtime_otimiz / 3, runtime_normal / 3 >> "tempo_triang.csv"
  }

  # Imprime cabeçalho para Banda de Memória
  print "#      MemoryBandwidth[MBytes/s]" > "mem_triang.csv"
  print "# N TriangularizaOtimiz Triangulariza" >> "mem_triang.csv"
  for (isize in graph) {
      if (!isize) continue

      printf "%s ", isize >> "mem_triang.csv"

      bandwidth_otimiz = graph[isize]["L3"]["Region TriangularizaOtimiz"]["L3 bandwidth [MBytes/s]"]
      bandwidth_normal = graph[isize]["L3"]["Region Triangulariza"]["L3 bandwidth [MBytes/s]"]

      printf "%f %f\n", bandwidth_otimiz, bandwidth_normal >> "mem_triang.csv"
  }

  # Imprime cabeçalho para Cache miss L2
  print "#      DataCacheMissRatio" > "cmiss_triang.csv"
  print "# N TriangularizaOtimiz Triangulariza" >> "cmiss_triang.csv"
  for (isize in graph) {
      if (!isize) continue

      printf "%s ", isize >> "cmiss_triang.csv"

      cmiss_otimiz = graph[isize]["L2CACHE"]["Region TriangularizaOtimiz"]["L2 miss ratio"]
      cmiss_normal = graph[isize]["L2CACHE"]["Region Triangulariza"]["L2 miss ratio"]

      printf "%f %f\n", cmiss_otimiz, cmiss_normal >> "cmiss_triang.csv"
  }

  # Imprime cabeçalho para Cache miss L2
  print "#      FLOPS_DP" > "flops_triang.csv"
  print "# N TriangularizaOtimiz Triangulariza" >> "flops_triang.csv"
  for (isize in graph) {
      if (!isize) continue

      printf "%s ", isize >> "flops_triang.csv"

      flops_otimiz = graph[isize]["FLOPS_DP"]["Region TriangularizaOtimiz"]["DP MFLOP/s"]
      flops_normal = graph[isize]["FLOPS_DP"]["Region Triangulariza"]["DP MFLOP/s"]

      printf "%f %f\n", flops_otimiz, flops_normal >> "flops_triang.csv"
  }
}
