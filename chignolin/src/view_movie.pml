# Visor INTERACTIVO de la trayectoria de Chignolin (estilo Folding@Home).
# Abre una ventana 3D y reproduce la simulación como película.
#
# Uso (en tu escritorio):
#     conda activate jic-folding
#     pymol src/view_movie.pml
#
# Para ver OTRA trayectoria, cambia el nombre del .dcd más abajo.

load results/prep/system.pdb, chig
load_traj results/trajectories/round0_000.dcd, chig, 1

remove not polymer      # quitar agua e iones (deja solo la proteína) -- robusto
remove hydrogens        # se ve más limpio

dss                     # asignar estructura secundaria (para el "cartoon")
hide everything
show cartoon
show sticks
color green, elem C
color red,   elem O
color blue,  elem N
set cartoon_transparency, 0.3
bg_color white

intra_fit chig          # alinear los frames (que la proteína no "vuele")
orient
set movie_fps, 15
# NOTA: no usamos 'mplay' (auto-reproducir desde script da un bug en PyMOL).
# Para ver la película: pulsa el botón  >  (Play), abajo a la derecha de la ventana.
