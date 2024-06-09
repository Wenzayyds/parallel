#include <iostream>
#include <mpi.h>
#include <sys/time.h>
#include <pmmintrin.h>
#include <omp.h>

using namespace std;

static const int N = 1000;
static const int thread_count = 4;

float arr[N][N];
float A[N][N];

void init_A(float arr[][N])
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            arr[i][j] = 0;
        }
        arr[i][i] = 1.0;
        for (int j = i + 1; j < N; j++)
            arr[i][j] = rand() % 100;
    }

    for (int i = 0; i < N; i++)
    {
        int k1 = rand() % N;
        int k2 = rand() % N;
        for (int j = 0; j < N; j++)
        {
            arr[i][j] += arr[0][j];
            arr[k1][j] += arr[k2][j];
        }
    }
}

void reset_A(float A[][N], float arr[][N])
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A[i][j] = arr[i][j];
}


void print_A(float A[][N])
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            cout << A[i][j] << " ";
        cout << endl;
    }
    cout << endl;
}



void LU(float A[][N], int rank, int num_proc)
{
    // 计算当前进程的前一个和后一个进程的编号
    int pre_proc = (rank + num_proc - 1) % num_proc;
    int next_proc = (rank + 1) % num_proc;

    // 遍历每一行进行LU分解
    for (int k = 0; k < N; k++)
    {
        // 检查当前行是否由本进程处理
        if (k % num_proc == rank)
        {
            // 对当前行进行规范化处理
            for (int j = k + 1; j < N; j++)
                A[k][j] /= A[k][k];
            A[k][k] = 1.0;

            // 将处理后的行发送给下一个进程
            MPI_Send(&A[k], N, MPI_FLOAT, next_proc, 2, MPI_COMM_WORLD);
        }
        else
        {
            // 接收前一个进程处理过的行
            MPI_Recv(&A[k], N, MPI_FLOAT, pre_proc, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // 如果当前行不是下一个进程要处理的行，继续转发给下一个进程
            if (k % num_proc != next_proc)
                MPI_Send(&A[k], N, MPI_FLOAT, next_proc, 2, MPI_COMM_WORLD);
        }

        // 更新当前进程负责的行
        for (int i = k + 1; i < N; i++)
        {
            if (i % num_proc == rank)
            {
                for (int j = k + 1; j < N; j++)
                    A[i][j] -= A[i][k] * A[k][j];
                A[i][k] = 0.0;
            }
        }
    }
}

void f_mpi()
{
    timeval t_start; // 定义开始时间结构体
    timeval t_end;   // 定义结束时间结构体

    int num_proc;    // 进程总数
    int rank;        // 当前进程的编号

    // 获取MPI环境中的进程总数和当前进程的编号
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // 如果是0号进程，执行特定的初始化和数据分发
    if (rank == 0)
    {
        reset_A(A, arr); // 重置矩阵A
        gettimeofday(&t_start, NULL); // 获取开始时间

        // 将矩阵的每一行根据行号模进程数分发给对应的进程
        for (int i = 0; i < N; i++)
        {
            int flag = i % num_proc;
            if (flag == rank)
                continue;
            else
                MPI_Send(&A[i], N, MPI_FLOAT, flag, 0, MPI_COMM_WORLD);
        }
        LU(A, rank, num_proc); // 执行LU分解

        // 接收其他进程处理完毕的数据
        for (int i = 0; i < N; i++)
        {
            int flag = i % num_proc;
            if (flag == rank)
                continue;
            else
                MPI_Recv(&A[i], N, MPI_FLOAT, flag, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        gettimeofday(&t_end, NULL); // 获取结束时间

        // 计算并输出总的时间消耗
        cout << "Pipeline MPI LU time cost: "
            << 1000 * (t_end.tv_sec - t_start.tv_sec) +
            0.001 * (t_end.tv_usec - t_start.tv_usec) << "ms" << endl;
    }
    else
    {
        // 非0号进程接收0号进程分发的数据
        for (int i = rank; i < N; i += num_proc)
        {
            MPI_Recv(&A[i], N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        LU(A, rank, num_proc); // 执行LU分解

        // 将处理完的数据发送回0号进程
        for (int i = rank; i < N; i += num_proc)
        {
            MPI_Send(&A[i], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        }
    }
}



void LU_opt(float A[][N], int rank, int num_proc)
{
    __m128 t1, t2, t3;
    int pre_proc = (rank + (num_proc - 1)) % num_proc;
    int next_proc = (rank + 1) % num_proc;
#pragma omp parallel num_threads(thread_count)
    for (int k = 0; k < N; k++)
    {
        if (int(k % num_proc) == rank)
        {
            float temp1[4] = { A[k][k], A[k][k], A[k][k], A[k][k] };
            t1 = _mm_loadu_ps(temp1);
            
#pragma omp for schedule(dynamic, 20)
            for (int j = k + 1; j < N - 3; j += 4)
            {
                t2 = _mm_loadu_ps(A[k] + j);
                t3 = _mm_div_ps(t2, t1);
                _mm_storeu_ps(A[k] + j, t3);
            }
            for (int j = N - N % 4; j < N; j++)
            {
                A[k][j] = A[k][j] / A[k][k];
            }
            A[k][k] = 1.0;
            MPI_Send(&A[k], N, MPI_FLOAT, next_proc, 2, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Recv(&A[k], N, MPI_FLOAT, pre_proc, 2,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (int(k % num_proc) != next_proc)
                MPI_Send(&A[k], N, MPI_FLOAT, next_proc, 2, MPI_COMM_WORLD);
        }
        for (int i = k + 1; i < N; i++)
        {
            if (int(i % num_proc) == rank)
            {
                float temp2[4] = { A[i][k], A[i][k], A[i][k], A[i][k] };
                t1 = _mm_loadu_ps(temp2);
#pragma omp for schedule(dynamic, 20)
                for (int j = k + 1; j <= N - 3; j += 4)
                {
                    t2 = _mm_loadu_ps(A[i] + j);
                    t3 = _mm_loadu_ps(A[k] + j);
                    t3 = _mm_mul_ps(t1, t3);
                    t2 = _mm_sub_ps(t2, t3);
                    _mm_storeu_ps(A[i] + j, t2);
                }
                for (int j = N - N % 4; j < N; j++)
                    A[i][j] = A[i][j] - A[i][k] * A[k][j];
                A[i][k] = 0;
            }
        }
    }
}

void f_mpi_opt()
{
    timeval t_start;
    timeval t_end;

    int num_proc;
    int rank;

    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        reset_A(A, arr);
        gettimeofday(&t_start, NULL);
        for (int i = 0; i < N; i++)
        {
            int flag = i % num_proc;
            if (flag == rank)
                continue;
            else
                MPI_Send(&A[i], N, MPI_FLOAT, flag, 0, MPI_COMM_WORLD);
        }
        LU_opt(A, rank, num_proc);
        for (int i = 0; i < N; i++)
        {
            int flag = i % num_proc;
            if (flag == rank)
                continue;
            else
                MPI_Recv(&A[i], N, MPI_FLOAT, flag, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        gettimeofday(&t_end, NULL);
        cout << "Pipeline MPI LU with SSE and OpenMP time cost: "
            << 1000 * (t_end.tv_sec - t_start.tv_sec) +
            0.001 * (t_end.tv_usec - t_start.tv_usec) << "ms" << endl;
        //print_A(A);
    }
    else
    {
        for (int i = rank; i < N; i += num_proc)
        {
            MPI_Recv(&A[i], N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        LU_opt(A, rank, num_proc);
        for (int i = rank; i < N; i += num_proc)
        {
            MPI_Send(&A[i], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        }
    }
}


int main()
{
    init_A(arr);

    MPI_Init(NULL, NULL);

    f_mpi();
    f_mpi_opt();
    MPI_Finalize();



}