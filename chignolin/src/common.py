"""Funciones compartidas por todos los scripts del pipeline.

Centraliza la lectura de config y la creación de objetos de OpenMM para que
`prepare_system.py`, `run_md.py` y `adaptive_sampling.py` usen exactamente la
misma configuración.
"""
import os
import yaml
import openmm as mm
from openmm import app, unit


def load_config(path="config.yaml"):
    """Carga config.yaml como un diccionario."""
    with open(path) as f:
        return yaml.safe_load(f)


def load_system(system_xml):
    """Lee un System de OpenMM previamente serializado."""
    with open(system_xml) as f:
        return mm.XmlSerializer.deserialize(f.read())


def make_integrator(cfg):
    """Crea el integrador de Langevin (termostato) según la config."""
    return mm.LangevinMiddleIntegrator(
        cfg["temperature_K"] * unit.kelvin,
        cfg["friction_per_ps"] / unit.picosecond,
        cfg["timestep_ps"] * unit.picoseconds,
    )


def make_simulation(topology, system, integrator, cfg):
    """Crea la Simulation en la plataforma elegida (CUDA por defecto)."""
    platform = mm.Platform.getPlatformByName(cfg.get("platform", "CUDA"))
    props = {}
    if cfg.get("platform", "CUDA") in ("CUDA", "OpenCL"):
        props["Precision"] = cfg.get("precision", "mixed")
    return app.Simulation(topology, system, integrator, platform, props)


def steps_for_ns(ns, cfg):
    """Convierte nanosegundos a número de pasos de integración."""
    return int(round(ns * 1000.0 / cfg["timestep_ps"]))


def frames_interval(cfg):
    """Cada cuántos pasos se guarda un frame."""
    return int(round(cfg["save_interval_ps"] / cfg["timestep_ps"]))
