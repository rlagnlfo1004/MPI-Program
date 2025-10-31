# MPI-program
분산 병렬 환경에서의 mpi 프로그램 작성

## 2D Heat Diffusion using MPI (MPI_Send/Recv Version)

https://github.com/rlagnlfo1004/mpi-program/blob/main/heat_diffusion/heat_diffusion.c

본 프로그램은 2차원 공간에서 열확산(heat diffusion) 과정을 시뮬레이션하는 MPI 프로그램입니다.

- 입력: 전역 격자 크기(NX, NY), 프로세스 격자(PX, PY), 반복 횟수(steps), 확산 계수(alpha)를 명령줄 인자로 받습니다.

- 출력: 최종 상태의 전체 온도 합(Global sum)과 최고 온도(Global max)를 계산하여 출력합니다.

- 주요 로직: 전체 격자를 PX * PY 개의 영역으로 분할(Domain Decomposition)하여 각 프로세스에 할당합니다. 각 프로세스는 MPI_Send와 MPI_Recv만을 사용하여 이웃 프로세스와 경계 데이터(Halo)를 교환하며, 5-point Stencil 계산을 통해 시간에 따른 열 분포를 시뮬레이션합니다.