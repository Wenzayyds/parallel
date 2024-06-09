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

void f_ordinary()
{
    reset_A(A, arr);
    timeval t_start;
    timeval t_end;
    gettimeofday(&t_start, NULL);

    for (int k = 0; k < N; k++)
    {
        for (int j = k + 1; j < N; j++)
        {
            A[k][j] = A[k][j] * 1.0 / A[k][k];
        }
        A[k][k] = 1.0;

        for (int i = k + 1; i < N; i++)
        {
            for (int j = k + 1; j < N; j++)
            {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
    gettimeofday(&t_end, NULL);
    cout << "ordinary time cost: "
        << 1000 * (t_end.tv_sec - t_start.tv_sec) +
        0.001 * (t_end.tv_usec - t_start.tv_usec) << "ms" << endl;
}



void LU(float A[][N], int rank, int num_proc)
{
    // 计算每个进程处理的块大小和剩余部分
    int block = N / num_proc;
    int remain = N % num_proc;

    // 计算当前进程处理的行的起始和结束索引
    int begin = rank * block;
    int end = (rank != num_proc - 1) ? (begin + block) : (begin + block + remain);
   
    // 对矩阵进行LU分解
    for (int k = 0; k < N; k++)
    {
        // 如果当前行由本进程处理
        if (k >= begin && k < end)
        {
            // 将主对角线元素下方的元素变为0，并将主对角线上的元素规范化为1
            for (int j = k + 1; j < N; j++)
                A[k][j] = A[k][j] / A[k][k];
            A[k][k] = 1.0;

            // 将当前行广播给其他进程
            for (int p = 0; p < num_proc; p++)
                if (p != rank)
                    MPI_Send(&A[k], N, MPI_FLOAT, p, 2, MPI_COMM_WORLD);

        }

        else
        {
            // 如果当前行不由本进程处理，则从负责该行的进程接收数据
            int cur_p = k / block;
             MPI_Recv(&A[k], N, MPI_FLOAT, cur_p, 2,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }

        // 更新当前进程负责的行
        for (int i = begin; i < end && i < N; i++)
        {
            if (i >= k + 1)
            {
                for (int j = k + 1; j < N; j++)
                    A[i][j] = A[i][j] - A[i][k] * A[k][j];
                A[i][k] = 0.0;
            }
        }
    }
}


void f_mpi()
{

    timeval t_start;
    timeval t_end;

    int num_proc;//������
    int rank;//ʶ����ý��̵�rank��ֵ��0~size-1

    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int block = N / num_proc;
    int remain = N % num_proc;

    //0�Ž��̡������񻮷�
    if (rank == 0)
    {
        reset_A(A, arr);
        gettimeofday(&t_start, NULL);
        //���񻮷�
        for (int i = 1; i < num_proc; i++)
        {
            if (i != num_proc - 1)
            {
                for (int j = 0; j < block; j++)
                    MPI_Send(&A[i * block + j], N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
            }
            else
            {
                for (int j = 0; j < block + remain; j++)
                    MPI_Send(&A[i * block + j], N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
            }
        }
        LU(A, rank, num_proc);
        //������0�Ž����Լ��������������������̴���֮��Ľ��
        for (int i = 1; i < num_proc; i++)
        {
            if (i != num_proc - 1)
            {
                for (int j = 0; j < block; j++)
                    MPI_Recv(&A[i * block + j], N, MPI_FLOAT, i, 1,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else
            {
                for (int j = 0; j < block + remain; j++)
                    MPI_Recv(&A[i * block + j], N, MPI_FLOAT, i, 1,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        gettimeofday(&t_end, NULL);
        cout << "Block MPI LU time cost: "
            << 1000 * (t_end.tv_sec - t_start.tv_sec) +
            0.001 * (t_end.tv_usec - t_start.tv_usec) << "ms" << endl;
        //print_A(A);
    }

    //��������
    else
    {
        //��0�Ž����Ƚ�������
        if (rank != num_proc - 1)
        {
            for (int j = 0; j < block; j++)
                MPI_Recv(&A[rank * block + j], N, MPI_FLOAT, 0, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else
        {
            for (int j = 0; j < block + remain; j++)
                MPI_Recv(&A[rank * block + j], N, MPI_FLOAT, 0, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        LU(A, rank, num_proc);
        //��0�Ž����������֮�󣬽�������ص�0�Ž���
        if (rank != num_proc - 1)
        {
            for (int j = 0; j < block; j++)
                MPI_Send(&A[rank * block + j], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        }
        else
        {
            for (int j = 0; j < block + remain; j++)
                MPI_Send(&A[rank * block + j], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        }
    }
}


void LU_opt(float A[][N], int rank, int num_proc)
{
    __m128 t1, t2, t3; // 定义SIMD向量变量
    int block = N / num_proc; // 计算每个进程处理的块大小
    int remain = N % num_proc; // 计算剩余的行数
    int begin = rank * block; // 计算当前进程处理的起始行
    int end = rank != num_proc - 1 ? begin + block : begin + block + remain; // 计算当前进程处理的结束行

    // 使用OpenMP并行化处理
#pragma omp parallel num_threads(thread_count),private(t1, t2, t3)
    for (int k = 0; k < N; k++)
    {
        if (k >= begin && k < end) // 当前进程负责处理的行
        {
            float temp1[4] = { A[k][k], A[k][k], A[k][k], A[k][k] }; // 将对角元素复制四次
            t1 = _mm_loadu_ps(temp1); // 加载对角元素到SIMD向量t1
#pragma omp for schedule(static)
            for (int j = k + 1; j < N - 3; j += 4) // 使用SIMD进行向量化计算
            {
                t2 = _mm_loadu_ps(A[k] + j);
                t3 = _mm_div_ps(t2, t1);
                _mm_storeu_ps(A[k] + j, t3);
            }
            for (int j = N - N % 4; j < N; j++) // 处理剩余的元素
            {
                A[k][j] = A[k][j] / A[k][k];
            }
            A[k][k] = 1.0; // 将对角元素设置为1
            for (int p = rank + 1; p < num_proc; p++)
                MPI_Send(&A[k], N, MPI_FLOAT, p, 2, MPI_COMM_WORLD); // 将处理后的行发送给其他进程
        }
        else
        {
            int cur_p = k / block; // 计算负责处理当前行的进程编号
            if (cur_p < rank)
                MPI_Recv(&A[k], N, MPI_FLOAT, cur_p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // 接收其他进程发送的行
        }
        for (int i = begin; i < end && i < N; i++) // 更新当前进程负责的行
        {
            if (i >= k + 1)
            {
                float temp2[4] = { A[i][k], A[i][k], A[i][k], A[i][k] }; // 将A[i][k]复制四次
                t1 = _mm_loadu_ps(temp2); // 加载A[i][k]到SIMD向量t1
#pragma omp for schedule(static)
                for (int j = k + 1; j <= N - 3; j += 4) // 使用SIMD进行向量化计算
                {
                    t2 = _mm_loadu_ps(A[i] + j); // 加载A[i][j]到SIMD向量t2
                    t3 = _mm_loadu_ps(A[k] + j); // 加载A[k][j]到SIMD向量t3
                    t3 = _mm_mul_ps(t1, t3); // 对t1中的每个元素乘以t3中的对应元素
                    t2 = _mm_sub_ps(t2, t3); // 对t2中的每个元素减去t3中的对应元素
                    _mm_storeu_ps(A[i] + j, t2); // 将计算结果存回A[i][j]
                }
                for (int j = N - N % 4; j < N; j++) // 处理剩余的元素
                    A[i][j] = A[i][j] - A[i][k] * A[k][j];
                A[i][k] = 0; // 将A[i][k]设置为0
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

    int block = N / num_proc;
    int remain = N % num_proc;

    if (rank == 0)
    {
        reset_A(A, arr);
        gettimeofday(&t_start, NULL);

        for (int i = 1; i < num_proc; i++)
        {
            if (i != num_proc - 1)
            {
                for (int j = 0; j < block; j++)
                    MPI_Send(&A[i * block + j], N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
            }
            else
            {
                for (int j = 0; j < block + remain; j++)
                    MPI_Send(&A[i * block + j], N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
            }
        }
        LU_opt(A, rank, num_proc);
        for (int i = 1; i < num_proc; i++)
        {
            if (i != num_proc - 1)
            {
                for (int j = 0; j < block; j++)
                    MPI_Recv(&A[i * block + j], N, MPI_FLOAT, i, 1,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else
            {
                for (int j = 0; j < block + remain; j++)
                    MPI_Recv(&A[i * block + j], N, MPI_FLOAT, i, 1,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        gettimeofday(&t_end, NULL);
        cout << "Block MPI LU with SSE and OpenMP time cost: "
            << 1000 * (t_end.tv_sec - t_start.tv_sec) +
            0.001 * (t_end.tv_usec - t_start.tv_usec) << "ms" << endl;
    }


    else
    {
        if (rank != num_proc - 1)
        {
            for (int j = 0; j < block; j++)
                MPI_Recv(&A[rank * block + j], N, MPI_FLOAT, 0, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else
        {
            for (int j = 0; j < block + remain; j++)
                MPI_Recv(&A[rank * block + j], N, MPI_FLOAT, 0, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        LU_opt(A, rank, num_proc);
        if (rank != num_proc - 1)
        {
            for (int j = 0; j < block; j++)
                MPI_Send(&A[rank * block + j], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        }
        else
        {
            for (int j = 0; j < block + remain; j++)
                MPI_Send(&A[rank * block + j], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        }
    }
}


int main()
{
    init_A(arr);

    //f_ordinary();

    MPI_Init(NULL, NULL);

    f_mpi();
    f_mpi_opt();
    MPI_Finalize();

}