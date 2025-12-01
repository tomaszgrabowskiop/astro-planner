#!/bin/bash
cd /Users/wdrodze/Public/astroplaner
source .venv/bin/activate
streamlit run RedCatPlanner.py --server.port 8501 --server.headless true