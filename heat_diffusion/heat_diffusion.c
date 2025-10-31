#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 입력 파라미터 처리 (NX, NY, PX, PY, steps, alpha)
    if (argc != 7) {
        if (rank == 0) printf("Usage: %s NX NY PX PY steps alpha\n", argv[0]);
        MPI_Finalize();
        return 1;
    }
    int NX = atoi(argv[1]);
    int NY = atoi(argv[2]);
    int PX = atoi(argv[3]);
    int PY = atoi(argv[4]);
    int steps = atoi(argv[5]);
    double alpha = atof(argv[6]);

    if (PX * PY != size) {
        if (rank == 0) printf("Error: PX * PY must be equal to the number of processes.\n");
        MPI_Finalize();
        return 1;
    }

    // 자신의 좌표 (my_px, my_py) 계산
    int my_px = rank % PX;
    int my_py = rank / PX;

    // 이웃 상, 하, 좌, 우 rank 계산 (경계 밖이면 MPI_PROC_NULL 사용)
    int north = (my_py > 0) ? (my_py - 1) * PX + my_px : MPI_PROC_NULL;
    int south = (my_py < PY - 1) ? (my_py + 1) * PX + my_px : MPI_PROC_NULL;
    int west = (my_px > 0) ? my_py * PX + (my_px - 1) : MPI_PROC_NULL;
    int east = (my_px < PX - 1) ? my_py * PX + (my_px + 1) : MPI_PROC_NULL;

    // Local Grid 크기 계산
    int base_nx = NX / PX;
    int rem_nx  = NX % PX;
    int local_nx = base_nx + (my_px < rem_nx ? 1 : 0);

    int base_ny = NY / PY;
    int rem_ny  = NY % PY;
    int local_ny = base_ny + (my_py < rem_ny ? 1 : 0);

    // Halo Cell을 위해 상하좌우 + 2 크기로 할당
    double *T_current = (double *)malloc((local_nx + 2) * (local_ny + 2) * sizeof(double));
    double *T_next = (double *)malloc((local_nx + 2) * (local_ny + 2) * sizeof(double));
    double **T_current_2d = malloc((local_ny + 2) * sizeof(double *));
    double **T_next_2d = (double **)malloc((local_ny + 2) * sizeof(double *));
    for (int i = 0; i < local_ny + 2; i++) {
        T_current_2d[i] = &T_current[i * (local_nx + 2)];
        T_next_2d[i] = &T_next[i * (local_nx + 2)];
    }

    // 모든 온도를 0으로 초기화
    for (int i = 0; i < (local_ny + 2) * (local_nx + 2); i++) {
        T_current[i] = 0.0;
        T_next[i] = 0.0;
    }

    // 열점(Hotspot) 설정 [0,0], 100도
    if (rank == 0) {
        T_current_2d[1][1] = 100.0;
    }

    // 좌/우 경계 통신을 위한 임시 버퍼
    double *send_buf = (double*)malloc(local_ny * sizeof(double));
    double *recv_buf = (double*)malloc(local_ny * sizeof(double));

    // A. 메인 steps 실행
    for (int t = 0; t < steps; t++) {
        // 전역 외곽 경계 조건 적용 (해당 경계 프로세스만)
        // Neumann (단열, zero heat-flux), 즉 바깥 halo 값은 인접 내부값을 복사
        if (north == MPI_PROC_NULL) { // 가장 상단 프로세스
            for (int j = 1; j <= local_nx; j++) T_current_2d[0][j] = T_current_2d[1][j];
        }
        if (south == MPI_PROC_NULL) { // 가장 하단 프로세스
            for (int j = 1; j <= local_nx; j++) T_current_2d[local_ny + 1][j] = T_current_2d[local_ny][j];
        }
        if (west == MPI_PROC_NULL) { // 가장 왼쪽 프로세스
            for (int i = 1; i <= local_ny; i++) T_current_2d[i][0] = T_current_2d[i][1];
        }
        if (east == MPI_PROC_NULL) { // 가장 오른쪽 프로세스
            for (int i = 1; i <= local_ny; i++) T_current_2d[i][local_nx + 1] = T_current_2d[i][local_nx];
        }

        // 경계값 교환(Halo exchange)
        MPI_Status status;

        // (1) 위로 보내고 아래로부터 받기
        MPI_Send(&T_current_2d[1][1], local_nx, MPI_DOUBLE, north, 0, MPI_COMM_WORLD);
        MPI_Recv(&T_current_2d[local_ny + 1][1], local_nx, MPI_DOUBLE, south, 0, MPI_COMM_WORLD, &status);

        // (2) 아래로 보내고 위로부터 받기
        MPI_Send(&T_current_2d[local_ny][1], local_nx, MPI_DOUBLE, south, 1, MPI_COMM_WORLD);
        MPI_Recv(&T_current_2d[0][1], local_nx, MPI_DOUBLE, north, 1, MPI_COMM_WORLD, &status);

        // (3) 왼쪽으로 보내고 오른쪽으로부터 받기 (버퍼 사용)
        for(int i = 0; i < local_ny; i++) send_buf[i] = T_current_2d[i+1][1]; // 버퍼로 데이터 복사
        MPI_Send(send_buf, local_ny, MPI_DOUBLE, west, 2, MPI_COMM_WORLD);
        MPI_Recv(recv_buf, local_ny, MPI_DOUBLE, east, 2, MPI_COMM_WORLD, &status);
        if (east != MPI_PROC_NULL) { // 이웃이 있을 때만
            for(int i = 0; i < local_ny; i++) T_current_2d[i+1][local_nx+1] = recv_buf[i];
        }

        // (4) 오른쪽으로 보내고 왼쪽으로부터 받기 (버퍼 사용)
        for(int i = 0; i < local_ny; i++) send_buf[i] = T_current_2d[i+1][local_nx]; // 버퍼로 데이터 복사
        MPI_Send(send_buf, local_ny, MPI_DOUBLE, east, 3, MPI_COMM_WORLD);
        MPI_Recv(recv_buf, local_ny, MPI_DOUBLE, west, 3, MPI_COMM_WORLD, &status);
        if (west != MPI_PROC_NULL) { // 이웃이 있을 때만
            for(int i = 0; i < local_ny; i++) T_current_2d[i+1][0] = recv_buf[i];
        }

        // 열확산 계산(5-point stencil)
        for (int i = 1; i <= local_ny; i++) {
            for (int j = 1; j <= local_nx; j++) {
                double N = T_current_2d[i - 1][j];
                double S = T_current_2d[i + 1][j];
                double W = T_current_2d[i][j - 1];
                double E = T_current_2d[i][j + 1];
                double C = T_current_2d[i][j];
                T_next_2d[i][j] = C + alpha * (N + S + W + E - 4 * C);
            }
        }

        // 더블 버퍼(Double Buffer)
        double *temp_buffer = T_current;
        T_current = T_next;
        T_next = temp_buffer;

        for (int i = 0; i < local_ny + 2; i++) {
            T_current_2d[i] = &T_current[i * (local_nx + 2)];
            T_next_2d[i] = &T_next[i * (local_nx + 2)];
        }
    }

    // B. 통계 계산
    double local_sum = 0.0;
    double local_max = 0.0;

    // 각 프로세스가 로컬 합/최댓값 계산
    for (int i = 1; i <= local_ny; i++) {
        for (int j = 1; j <= local_nx; j++) {
            local_sum += T_current_2d[i][j];
            if (T_current_2d[i][j] > local_max) {
                local_max = T_current_2d[i][j];
            }
        }
    }

    if (rank == 0) { // rank 0이 모든 결과를 수집
        double global_sum = local_sum;
        double global_max = local_max;

        MPI_Status status;
        for (int i = 1; i < size; ++i) {
            MPI_Recv(&local_sum, 1, MPI_DOUBLE, i, 10, MPI_COMM_WORLD, &status);
            MPI_Recv(&local_max, 1, MPI_DOUBLE, i, 11, MPI_COMM_WORLD, &status);
            global_sum += local_sum;
            if (local_max > global_max) {
                global_max = local_max;
            }
        }
        // 최종 결과 출력
        printf("Global sum: %f\n", global_sum);
        printf("Global max: %f\n", global_max);
    } else {
        // rank 0에게 자신의 local_sum과 local_max 전송
        MPI_Send(&local_sum, 1, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
        MPI_Send(&local_max, 1, MPI_DOUBLE, 0, 11, MPI_COMM_WORLD);
    }

    free(T_current);
    free(T_next);
    free(T_current_2d);
    free(T_next_2d);
    free(send_buf);
    free(recv_buf);

    MPI_Finalize();
    return 0;
}