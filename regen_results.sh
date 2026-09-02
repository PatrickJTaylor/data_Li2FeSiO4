#!/bin/bash

# SPDX-FileCopyrightText: Copyright (C) 2026 Patrick J. Taylor
# SPDX-License-Identifier: GPL-3.0-or-later

set -o errexit
set -o nounset
set -o pipefail

CONTAINER_ID="regen_results_$(date +%s)_$RANDOM"

log() {
  printf "%b\n" "$*"
}

silence_until_error() {
  local out

  out=$("${@}" 2>&1) || {
    log "${out}" >&2
    exit 1
  }
}

contains() {
  local target_value=$1
  shift
  local valid_values=("${@}")
  for valid_value in "${valid_values[@]}"; do
    [[ "${valid_value}" == "${target_value}" ]] && return 0
  done

  return 1
}

spinner() {
  local text=$1
  local pid=$2
  local spinner_string='|/-\'

  local i=0
  while kill -0 "${pid}" 2>/dev/null; do
    i=$(( (i+1) %4 ))
    printf "\r%s %c " "${text}" "${spinner_string:$i:1}" >&2
    sleep 0.1
  done

  wait "${pid}" 2>/dev/null

  printf "\r%s    \n" "${text}" >&2
}

cleanup_container() {
  log "Cleaning up Docker container ${CONTAINER_ID}"

  docker rm --force "${CONTAINER_ID}" >/dev/null 2>&1
}

print_help() {
  cat << E0F
Usage: $0 [options] <all | figures | other_results | main | supp | <result_id> ...>

Regenerate the figures/results for the manuscript:
"Multi-Step Oxygen Redox Mechanism in the Polyanionic Lithium-Rich Cathode Li2FeSiO4"

Arguments:
  all                Regenerate everything (all figures + discussed-in-text results)
  figures            Regenerate all figures (main text + supplementary material, Fig. 2-S3)
  other_results      Regenerate all discussed-in-text results (oxi_states, thermo, fe_coord etc.)
  main               Regenerate only figures from the main text (Fig. 2-7)
  supp               Regenerate only figures from the supplementary material (Fig. S1-S3)
  <result_id> ...    Regenerate a specific set of results e.g., fig_2, fig_S1 and thermo

Possible values for <result_id> are:
  Figures:
    - fig_2         (Fe and O PDOS)
    - fig_3         (magmoms + ICOBIs at x = 2, 1, 0)
    - fig_5         (AIMD total energy evolution)
    - fig_6         (RDFs for selected MD frames)
    - fig_7         (oxidation states for selected MD frames)
    - fig_S1        (Li_xFeSiO4 polymorph stability at x = 2, 0)
    - fig_S2        (convex hull and voltage curve)
    - fig_S3        (Fe magmom and mean Fe-O ICOBI distributions by oxidation state)

  Other results (discussed in the text but not directly associated with a figure):
    - oxi_states    (oxidation states at x = 2, 1, 0)
    - thermo        (thermodynamic stability of FeSiO4)
    - fe_coord      (Fe coordination environments by oxidation state, selected MD frames)
    - fe_fe_pairs   (short Fe-Fe distances vs coordination and oxidation state)
    - o2_tracing    (O-O dimers and the Fe coordination of dimer-forming O atoms)
    - dimer_o_coord (full cation coordination of the 500 K peroxide-forming O atoms)

Options:
  -h, --help	Print this help message and exit

Examples:
  $0 all
  $0 figures
  $0 fig_7
  $0 fig_2 fig_3 fig_S1 fig_S2
  $0 fig_2 fig_5 oxi_states thermo
E0F
}

RESULT_IDS=("fig_2" "fig_3" "fig_5" "fig_6" "fig_7" "fig_S1" "fig_S2" "fig_S3" "oxi_states" "thermo" "fe_coord" "fe_fe_pairs" "o2_tracing" "dimer_o_coord")
FIG_IDS=("${RESULT_IDS[@]:0:8}")
TEXT_IDS=("${RESULT_IDS[@]:8}")
MAIN_IDS=("${FIG_IDS[@]:0:5}")
SUPP_IDS=("${FIG_IDS[@]:5}")
# Results that depend on Wannier-assigned oxidation states for the selected MD frames
NEEDS_MD_OX_STATES=("fig_7" "fig_S3" "fe_coord" "fe_fe_pairs" "o2_tracing" "dimer_o_coord")

if [[ $# -eq 0 ]] || [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
  print_help
  exit 0
fi

command -v docker >/dev/null 2>&1 || {
  log 'Error: docker not found on PATH' >&2
  exit 1
}

if [[ "$1" == "all" ]]; then
  RESULTS=("${RESULT_IDS[@]}")
elif [[ "$1" == "figures" ]]; then
  RESULTS=("${FIG_IDS[@]}")
elif [[ "$1" == "other_results" ]]; then
  RESULTS=("${TEXT_IDS[@]}")
elif [[ "$1" == "main" ]]; then
  RESULTS=("${MAIN_IDS[@]}")
elif [[ "$1" == "supp" ]]; then
  RESULTS=("${SUPP_IDS[@]}")
else
  RESULTS=("$@")

  for result_id in "${RESULTS[@]}"; do
    if ! contains "${result_id}" "${RESULT_IDS[@]}"; then
      log "Error: ${result_id} is not a valid result identifier" >&2
      exit 1
    fi
  done
fi

log 'Regenerating results:' "${RESULTS[@]}" "\n"

log 'Building Docker image with --tag result_regeneration'
silence_until_error docker build --tag result_regeneration .

log "Initialising Docker container with --name ${CONTAINER_ID}\n"
docker run --detach --name "${CONTAINER_ID}" result_regeneration sleep infinity >/dev/null

trap cleanup_container EXIT INT TERM

for result_id in "${RESULTS[@]}"; do
  if contains "${result_id}" "${NEEDS_MD_OX_STATES[@]}"; then
    docker exec "${CONTAINER_ID}" python "scripts/assign_md_oxidation_states.py" &
    spinner "Assigning Wannier oxidation states for selected MD frames (prerequisite of: ${result_id})" "${!}"

    printf "\n"
    break
  fi
done

for result_id in "${RESULTS[@]}"; do
  if [[ "${result_id}" == "oxi_states" ]]; then
    printf "\n"
    log 'Assigning Wannier oxidation states for x = 2, 1, 0:'
    docker exec "${CONTAINER_ID}" python "scripts/assign_delithiation_oxidation_states.py"
  elif [[ "${result_id}" == "thermo" ]]; then
    printf "\n"
    log 'Computing the thermodynamic stability of FeSiO4:'
    docker exec "${CONTAINER_ID}" python "scripts/thermodynamic_stability.py"
  elif [[ "${result_id}" == "fe_coord" ]]; then
    printf "\n"
    log 'Analysing Fe coordination environments by oxidation state for selected MD frames:'
    docker exec "${CONTAINER_ID}" python "scripts/fe_coordination_by_oxidation_state.py"
  elif [[ "${result_id}" == "fe_fe_pairs" ]]; then
    printf "\n"
    log 'Analysing short Fe-Fe distances for selected MD frames:'
    docker exec "${CONTAINER_ID}" python "scripts/short_fe_fe_distances.py"
  elif [[ "${result_id}" == "o2_tracing" ]]; then
    printf "\n"
    log 'Tracing O-O dimer formation through selected MD frames:'
    docker exec "${CONTAINER_ID}" python "scripts/trace_o2_forming_oxygens.py"
  elif [[ "${result_id}" == "dimer_o_coord" ]]; then
    printf "\n"
    log 'Analysing the cation coordination of the 500 K peroxide-forming O atoms:'
    docker exec "${CONTAINER_ID}" python "scripts/dimer_oxygen_coordination.py"
  else
    fig_num="${result_id:4}"

    docker exec "${CONTAINER_ID}" python "scripts/figure_${fig_num}.py" &
    spinner "Regenerating Fig. ${fig_num}" "${!}"
  fi
done

printf "\n"

log "Copying regenerated figures from ${CONTAINER_ID} -> ./figures/regenerated"
silence_until_error docker cp "${CONTAINER_ID}":/work/figures/regenerated ./figures
