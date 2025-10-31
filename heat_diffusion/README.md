# 2D Heat Diffusion using MPI (MPI_Send/Recv Version)

<br>

## 1. 프로그램 개요

본 프로그램은 2차원 공간에서 열확산(heat diffusion) 과정을 시뮬레이션하는 MPI 프로그램입니다.

- 입력: 전역 격자 크기(NX, NY), 프로세스 격자(PX, PY), 반복 횟수(steps), 확산 계수(alpha)를 명령줄 인자로 받습니다.

- 출력: 최종 상태의 전체 온도 합(Global sum)과 최고 온도(Global max)를 계산하여 출력합니다.

- 주요 로직: 전체 격자를 PX * PY 개의 영역으로 분할(Domain Decomposition)하여 각 프로세스에 할당합니다. 각 프로세스는 MPI_Send와 MPI_Recv만을 사용하여 이웃 프로세스와 경계 데이터(Halo)를 교환하며, 5-point Stencil 계산을 통해 시간에 따른 열 분포를 시뮬레이션합니다.

### 1.1. 제약 조건

- MPI_Send / MPI_Recv Point-to-Point 통신만 사용합니다.

- Collective Communication (MPI_Bcast, MPI_Reduce 등), MPI_Sendrecv, Non-blocking 통신 (MPI_Isend, MPI_Irecv)은 사용하지 않습니다.

- Single Machine에서 실행을 가정합니다.

<br/>

## 2. Build 및 Run Command

### 2.1. Build Command

```
mpicc -o heat_diffusion heat_diffusion.c
```
- `mpicc`: MPI 프로그램을 컴파일하기 위한 C 컴파일러 래퍼(wrapper)

### 2.2. Run Command
```
mpirun -np 16 ./hw1 512 512 4 4 100 0.23
```
- `mpirun`: MPI 프로그램을 실행하는 명령어
- `-np 16`: 사용할 총 프로세스의 개수 (PX * PY와 일치해야 함)
- `512 512 4 4 100 0.23`: 6개의 인자 (NX, NY, PX, PY, steps, alpha)


<br>

## 3. 코드 설명 (핵심 로직)

### 3.1. 초기 설정 및 Decomposition

- 입력 파라미터 처리: main 함수의 argv를 통해 6개의 인자(NX, NY, PX, PY, steps, alpha)를 받습니다. Rank 0이 인자 개수와 PX * PY == size를 검증합니다.

- Decomposition 및 로컬 그리드 설정: 전체 격자를 PX * PY개의 사각형 영역으로 분할합니다. 각 프로세스는 rank를 이용해 2차원 좌표 (my_px, my_py)를 계산합니다.

- Load Balancing: NX % PX 또는 NY % PY의 나머지를 앞 순위의 프로세스들이 하나씩 더 담당하도록 local_nx, local_ny를 결정합니다.

- 이웃 프로세스 식별: MPI_PROC_NULL을 활용합니다. 각 프로세스는 상/하/좌/우 이웃의 rank를 계산하며, 경계에 위치하여 이웃이 없는 경우 해당 rank를 MPI_PROC_NULL로 설정하여 if 조건문 없이 통신 로직을 단순화합니다.

### 3.2. 메모리 할당 및 초기화

- 더블 버퍼링 (Double Buffering): 현재 스텝(T_current)과 다음 스텝(T_next)을 위한 두 개의 독립된 버퍼를 할당합니다.

- Halo Cell: 각 로컬 격자 크기에 상하좌우 1칸씩의 Halo Cell(Ghost Cell)을 추가하여 (local_ny + 2) x (local_nx + 2) 크기의 배열을 동적으로 할당합니다.

- 2D 포인터 매핑: 계산 편의성을 위해 [row][col] 형태로 접근할 수 있도록 2차원 포인터 배열(T_current_2d, T_next_2d)로 매핑합니다.

- 초기 조건: 모든 셀을 0.0으로 초기화 후, Rank 0의 로컬 [1][1]에 100.0의 열점(Hotspot)을 설정합니다.

- 통신 버퍼: 좌/우 경계의 비연속적(non-contiguous) 데이터 전송을 위한 1차원 임시 버퍼(send_buf, recv_buf)를 local_ny 크기로 할당합니다.

### 3.3. 메인 steps 실행

1. 전역 외곽 경계 조건 적용: (해당 경계 프로세스만) 전역 격자의 가장자리에 위치한 프로세스들은 자신의 외부 halo 셀 값을 인접한 내부 셀 값으로 복제합니다.

2. 경계값 교환 (Halo Exchange):
  - 교착 상태(Deadlock) 방지: 순환 대기(Circular Wait)를 피하기 위해 송/수신 방향을 교차시키는 방식(e.g., (1)상(send)/하(recv) -> (2)하(send)/상(recv) -> (3)좌(send)/우(recv) -> (4)우(send)/좌(recv))을 사용합니다.
  - 수직 (상/하) 교환: 연속적인 메모리(행)이므로 직접 MPI_Send/MPI_Recv를 수행합니다.
  - 수평 (좌/우) 교환: 비연속적인 메모리(열)이므로, 임시 버퍼(send_buf, recv_buf)를 사용합니다. send_buf로 데이터를 복사(Packing)하여 송신하고, recv_buf로 수신한 뒤 다시 Halo 영역에 복사(Unpacking)합니다.

3. 열확산 계산 (5-point Stencil): T_current_2d를 읽어 5-point stencil 계산을 수행하고, 결과를 T_next_2d에 저장합니다.

4. 더블 버퍼링: T_current와 T_next의 포인터를 교체(swap)하고, 2차원 포인터 배열을 재설정합니다.

### 3.4. 통계 계산 (P2P 방식)

- 로컬 계산: 모든 steps 종료 후, 각 프로세스가 T_current 버퍼를 기준으로 local_sum과 local_max를 계산합니다.

- 데이터 취합:
  - rank > 0인 모든 프로세스가 local_sum과 local_max를 Rank 0에게 MPI_Send로 전송합니다.
  - rank 0은 for 루프를 돌며 rank 1부터 size-1까지의 모든 프로세스로부터 MPI_Recv를 호출하여 값을 수신합니다.

- 최종 계산 및 출력: Rank 0이 수신한 값들을 자신의 로컬 값에 누적/비교하여 최종 global_sum과 global_max를 계산하고 출력합니다.

### 3.5. 자원 해제

free()를 통해 동적으로 할당된 모든 메모리(데이터 버퍼, 2D 포인터 배열, 통신 버퍼)를 해제하고 MPI_Finalize()를 호출하여 MPI 환경을 종료합니다.

<br>

## 4. 코드 개선

### 4.1. MPI_PROC_NULL의 존재
- 초기 문제: 경계 프로세스(e.g., my_px == 0)를 처리하기 위해 복잡한 if-else 분기문을 사용했습니다.

- 경계 처리 로직 단순화 방법으로 MPI_PROC_NULL 상수를 사용했습니다.
  - 적용: MPI_PROC_NULL을 이웃 rank로 사용하면, 해당 MPI_Send/MPI_Recv 호출이 즉시 반환(no-op)됩니다. 이로 인해 모든 프로세스가 경계 여부와 관계없이 동일한 통신 로직을 사용하게 되어 코드가 매우 간결해졌습니다.

### 4.2. 비연속적 데이터 (좌/우 경계) 전송 방법

- 초기 문제: C언어 2D 배열은 행(row) 우선 저장되어, 열(column) 데이터가 메모리상에 비연속적입니다. 상/하 경계는 연속적이나 좌/우 경계 전송이 비효율적이었습니다.

- 1) MPI 파생 데이터 타입, 2) 임시 버퍼(temporary buffer) 사용의 두 가지 방법을 통해 이를 개선하였습니다.
  - 적용: 과제 요구사항(MPI_Send/MPI_Recv만 사용)을 고려하여 임시 버퍼 방식을 채택했습니다. (Packing) send_buf에 열 데이터를 복사하여 전송하고, (Unpacking) recv_buf로 수신 후 다시 halo 셀에 복사하는 로직을 직접 구현했습니다.

<br>

## 5. 참고한 사이트

Generative AI를 통해 MPI_PROC_NULL의 존재를 알게 된 후, MPI 표준 문서와 코드 예시들을 통해 정확한 동작 원리를 참고하였습니다.

- MPI-1.1 Standard
  - MPI_PROC_NULL이 통신에 사용될 때의 동작(즉시 성공적 완료)을 정의하는 공식 MPI 표준 문서.
  - 링크: https://www.mpi-forum.org/docs/mpi-1.1/mpi-11-html/node53.html

- RookieHPC
  - MPI_PROC_NULL의 개념과 사용 예제를 확인한 사이트. 경계 조건 처리를 간결하게 하는 예제 포함.
  - 링크: https://rookiehpc.org/mpi/docs/mpi_proc_null/index.html

<br>

## 6. 실험 결과

- `alpha` (확산 계수): 0.23
- `process` 수: 16
- `PX, PY`: 4

- | 격자 크기 | Step | 결과 (Global sum) | 결과 (Global max) |
  | :--- | :--- | :--- | :--- |
  | 512 X 512 | 100 | 100.000000 | 1.362933 |
  | 1024 X 1024 | 200 | 100.000000 | 0.686680 |
  | 2048 X 2048 | 400 | 100.000000 | 0.344659 |
  | 4096 X 4096 | 1000 | 100.000000 | 0.138182 |