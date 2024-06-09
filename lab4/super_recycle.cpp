# include <sys/time.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <mpi.h>
#include <pmmintrin.h>
#include <omp.h>
using namespace std;


static const int thread_count = 4;

/*
unsigned int Act[8399][264] = { 0 };
unsigned int Pas[8399][264] = { 0 };

const int Num = 263;
const int pasNum = 4535;
const int lieNum = 8399;
*/



//
//unsigned int Act[23045][722] = { 0 };
//unsigned int Pas[23045][722] = { 0 };
//
//const int Num = 721;
//const int pasNum = 14325;
//const int lieNum = 23045;


//
//unsigned int Act[37960][1188] = { 0 };
//unsigned int Pas[37960][1188] = { 0 };
//
//const int Num = 1187;
//const int pasNum = 14921;
//const int lieNum = 37960;



unsigned int Act[43577][1363] = { 0 };
unsigned int Pas[54274][1363] = { 0 };

const int Num = 1362;
const int pasNum = 54274;
const int lieNum = 43577;




void init_A()
{
    unsigned int a;
    ifstream infile("act3.txt");
    char fin[10000] = { 0 };
    int index;
    while (infile.getline(fin, sizeof(fin)))
    {
        std::stringstream line(fin);
        int biaoji = 0;
        while (line >> a)
        {
            if (biaoji == 0)
            {
                index = a;
                biaoji = 1;
            }
            int k = a % 32;
            int j = a / 32;

            int temp = 1 << k;
            Act[index][Num - 1 - j] += temp;
            Act[index][Num] = 1;
        }
    }
}

void init_P()
{
    unsigned int a;
    ifstream infile("pas3.txt");
    char fin[10000] = { 0 };
    int index = 0;
    while (infile.getline(fin, sizeof(fin)))
    {
        std::stringstream line(fin);
        int biaoji = 0;
        while (line >> a)
        {
            if (biaoji == 0)
            {
                Pas[index][Num] = a;
                biaoji = 1;
            }

            int k = a % 32;
            int j = a / 32;

            int temp = 1 << k;
            Pas[index][Num - 1 - j] += temp;
        }
        index++;
    }
}



void f_ordinary()
{
    timeval t_start;
    timeval t_end;
    gettimeofday(&t_start, NULL);

    bool sign;
    do
    {
        //---------------------------begin-------------------------------------

        int i;
        for (i = lieNum - 1; i - 8 >= -1; i -= 8)
        {
            for (int j = 0; j < pasNum; j++)
            {
                while (Pas[j][Num] <= i && Pas[j][Num] >= i - 7)
                {
                    int index = Pas[j][Num];

                    if (Act[index][Num] == 1)
                    {
                        for (int k = 0; k < Num; k ++)
                        {
                            Pas[j][k] = Pas[j][k] ^ Act[index][k];
                        }

                        int num = 0, S_num = 0;
                        for (num = 0; num < Num; num++)
                        {
                            if (Pas[j][num] != 0)
                            {
                                unsigned int temp = Pas[j][num];
                                while (temp != 0)
                                {
                                    temp = temp >> 1;
                                    S_num++;
                                }
                                S_num += num * 32;
                                break;
                            }
                        }
                        Pas[j][Num] = S_num - 1;
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }

        for (i = i + 8; i >= 0; i--)
        {

            for (int j = 0; j < pasNum; j++)
            {
                while (Pas[j][Num] == i)
                {
                    if (Act[i][Num] == 1)
                    {
                        for (int k = 0; k < Num; k ++)
                        {
                            Pas[j][k] = Pas[j][k] ^ Act[i][k];
                        }

                        int num = 0, S_num = 0;
                        for (num = 0; num < Num; num++)
                        {
                            if (Pas[j][num] != 0)
                            {
                                unsigned int temp = Pas[j][num];
                                while (temp != 0)
                                {
                                    temp = temp >> 1;
                                    S_num++;
                                }
                                S_num += num * 32;
                                break;
                            }
                        }
                        Pas[j][Num] = S_num - 1;

                    }
                    else
                    {
                        break;
                    }
                }
            }
        }

        //----------------------------------end--------------------------------



        //检查元素，然后判断是否继续
        sign = false;
        for (int i = 0; i < pasNum; i++)
        {
            //找到第i个非零元素的位置
            int temp = Pas[i][Num];
            if (temp == -1)
            {
                //说明已经变为零元素了
                continue;
            }

            //判断对应的行元素是否为空，不为空，则继续
            if (Act[temp][Num] == 0)
            {
                //变为行元素
                for (int k = 0; k < Num; k++)
                    Act[temp][k] = Pas[i][k];
                //行元素变为零
                Pas[i][Num] = -1;
                //标志bool为true，说明还需要继续
                sign = true;
            }
        }

    } while (sign == true);


    gettimeofday(&t_end, NULL);
    cout << "ordinary time cost: "
        << 1000 * (t_end.tv_sec - t_start.tv_sec) +
        0.001 * (t_end.tv_usec - t_start.tv_usec) << "ms" << endl;
}


void f_ordinary1()
{
    timeval t_start;
    timeval t_end;
    gettimeofday(&t_start, NULL);

    int i;
    for (i = lieNum - 1; i - 8 >= -1; i -= 8)
    {
        //每次处理8个行元素，范围是从 i-7 到 i

        for (int j = 0; j < pasNum; j++)
        {
            //对4535个行元素，没有在这个范围内的
            while (Pas[j][Num] <= i && Pas[j][Num] >= i - 7)
            {
                int index = Pas[j][Num];
                if (Act[index][Num] == 1)//行元素不为空
                {
                    //Pas[j][]和Act[在Pas[j][18]中][]异或
                    for (int k = 0; k < Num; k++)
                    {
                        Pas[j][k] = Pas[j][k] ^ Act[index][k];
                    }

                    //计算Pas[j][18]的最小值
                    //计算之后的结果存储在Pas[j][18]中，用于下次while循环
                    //计算之后Pas[j][ ]变为零
                    int num = 0, S_num = 0;
                    for (num = 0; num < Num; num++)
                    {
                        if (Pas[j][num] != 0)
                        {
                            unsigned int temp = Pas[j][num];
                            while (temp != 0)
                            {
                                temp = temp >> 1;
                                S_num++;
                            }
                            S_num += num * 32;
                            break;
                        }
                    }
                    Pas[j][Num] = S_num - 1;

                }
                else//行元素为空
                {
                    //Pas[j][]变为行元素
                    for (int k = 0; k < Num; k++)
                        Act[index][k] = Pas[j][k];

                    Act[index][Num] = 1;//变为行元素非空
                    break;
                }

            }
        }
    }


    for (int i = lieNum % 8 - 1; i >= 0; i--)
    {
        //每次处理1个行元素，范围是单个的i

        for (int j = 0; j < pasNum; j++)
        {
            //对53个行元素，没有在这个范围内的
            while (Pas[j][Num] == i)
            {
                if (Act[i][Num] == 1)//行元素不为空
                {
                    //Pas[j][]和Act[i][]异或
                    for (int k = 0; k < Num; k++)
                    {
                        Pas[j][k] = Pas[j][k] ^ Act[i][k];
                    }

                    //计算Pas[j][18]的最小值
                    //计算之后的结果存储在Pas[j][18]中，用于下次while循环
                    //计算之后Pas[j][ ]变为零
                    int num = 0, S_num = 0;
                    for (num = 0; num < Num; num++)
                    {
                        if (Pas[j][num] != 0)
                        {
                            unsigned int temp = Pas[j][num];
                            while (temp != 0)
                            {
                                temp = temp >> 1;
                                S_num++;
                            }
                            S_num += num * 32;
                            break;
                        }
                    }
                    Pas[j][Num] = S_num - 1;

                }
                else//行元素为空
                {
                    //Pas[j][]变为行元素
                    for (int k = 0; k < Num; k++)
                        Act[i][k] = Pas[j][k];

                    Act[i][Num] = 1;//变为行元素非空
                    break;
                }
            }
        }
    }


    gettimeofday(&t_end, NULL);
    cout << "ordinary time cost: "
        << 1000 * (t_end.tv_sec - t_start.tv_sec) +
        0.001 * (t_end.tv_usec - t_start.tv_usec) << "ms" << endl;
}


void super(int rank, int num_proc)
{
    //并行处理矩阵行元素------------------------------------------------------
       //---------------------------begin-------------------------------------


    int i;
    
    //每次处理8个行元素，范围是从 i-7 到 i
#pragma omp parallel num_threads(thread_count) 
    for (i = lieNum - 1; i - 8 >= -1; i -= 8)
    {
#pragma omp for schedule(dynamic,20)     
        for (int j = 0; j < pasNum; j++)
        {
            //当前处理自己分配的行元素。其他的忽略
            if (int(j % num_proc) == rank)
            {

                while (Pas[j][Num] <= i && Pas[j][Num] >= i - 7)
                {
                    int index = Pas[j][Num];

                    if (Act[index][Num] == 1)//行元素不为空
                    {
                        //****************************SIMD******************************
                        ////Pas[j][]和Act[在Pas[j][x]中][]异或
                        //for (int k = 0; k < Num; k++)
                        //{
                        //    Pas[j][k] = Pas[j][k] ^ Act[index][k];
                        //}
                        int k;
                        __m128 va_Pas, va_Act;
                        for (k = 0; k + 4 <= Num; k += 4)
                        {
                            //Pas[j][k] = Pas[j][k] ^ Act[index][k];
                            va_Pas = _mm_loadu_ps((float*)&(Pas[j][k]));
                            va_Act = _mm_loadu_ps((float*)&(Act[index][k]));

                            va_Pas = _mm_xor_ps(va_Pas, va_Act);
                            _mm_store_ss((float*)&(Pas[j][k]), va_Pas);
                        }

                        for (; k < Num; k++)
                        {
                            Pas[j][k] = Pas[j][k] ^ Act[index][k];
                        }
                        //***************************SIMD******************************


                        //计算Pas[j][18]的最小值
                        //计算之后的结果存储在Pas[j][18]中，用于下次while循环
                        //计算之后Pas[j][ ]变为零
                        int num = 0, S_num = 0;
                        for (num = 0; num < Num; num++)
                        {
                            if (Pas[j][num] != 0)
                            {
                                unsigned int temp = Pas[j][num];
                                while (temp != 0)
                                {
                                    temp = temp >> 1;
                                    S_num++;
                                }
                                S_num += num * 32;
                                break;
                            }
                        }
                        Pas[j][Num] = S_num - 1;
                    }
                    else//行元素为空
                    {
                        break;
                    }
                }
            }
        }
    }

#pragma omp parallel num_threads(thread_count) 
    for (int i = lieNum % 8 - 1; i >= 0; i--)
    {
        //每次处理1个行元素，范围是单个的i
#pragma omp for schedule(dynamic,20)
        for (int j = 0; j < pasNum; j++)
        {
            //当前处理自己分配的行元素。其他的忽略
            if (int(j % num_proc) == rank)
            {
                while (Pas[j][Num] == i)
                {
                    if (Act[i][Num] == 1)//行元素不为空
                    {
                         //****************************SIMD******************************
                        ////Pas[j][]和Act[在Pas[j][x]中][]异或
                        //for (int k = 0; k < Num; k++)
                        //{
                        //    Pas[j][k] = Pas[j][k] ^ Act[i][k];
                        //}
                        int k;
                        __m128 va_Pas, va_Act;
                        for (k = 0; k + 4 <= Num; k += 4)
                        {
                            //Pas[j][k] = Pas[j][k] ^ Act[index][k];
                            va_Pas = _mm_loadu_ps((float*)&(Pas[j][k]));
                            va_Act = _mm_loadu_ps((float*)&(Act[i][k]));

                            va_Pas = _mm_xor_ps(va_Pas, va_Act);
                            _mm_store_ss((float*)&(Pas[j][k]), va_Pas);
                        }

                        for (; k < Num; k++)
                        {
                            Pas[j][k] = Pas[j][k] ^ Act[i][k];
                        }
                        //***************************SIMD******************************


                        //计算Pas[j][18]的最小值
                        //计算之后的结果存储在Pas[j][18]中，用于下次while循环
                        //计算之后Pas[j][ ]变为零
                        int num = 0, S_num = 0;
                        for (num = 0; num < Num; num++)
                        {
                            if (Pas[j][num] != 0)
                            {
                                unsigned int temp = Pas[j][num];
                                while (temp != 0)
                                {
                                    temp = temp >> 1;
                                    S_num++;
                                }
                                S_num += num * 32;
                                break;
                            }
                        }
                        Pas[j][Num] = S_num - 1;

                    }
                    else//行元素为空
                    {
                        break;
                    }
                }
            }
        }
    }
    //----------------------------------end--------------------------------

}


void f_mpi()
{

    int num_proc;
    int rank;

    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        timeval t_start;
        timeval t_end;

        gettimeofday(&t_start, NULL);
        int sign;
        do
        {
            for (int i = 0; i < pasNum; i++)
            {
                int flag = i % num_proc;
                if (flag == rank)
                    continue;
                else
                    MPI_Send(&Pas[i], Num + 1, MPI_FLOAT, flag, 0, MPI_COMM_WORLD);
            }
            super(rank, num_proc);
            for (int i = 0; i < pasNum; i++)
            {
                int flag = i % num_proc;
                if (flag == rank)
                    continue;
                else
                    MPI_Recv(&Pas[i], Num + 1, MPI_FLOAT, flag, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            sign = 0;
            for (int i = 0; i < pasNum; i++)
            {
                int temp = Pas[i][Num];
                if (temp == -1)
                {
                    continue;
                }
                if (Act[temp][Num] == 0)
                {
                    for (int k = 0; k < Num; k++)
                        Act[temp][k] = Pas[i][k];
                    Pas[i][Num] = -1;
                    sign = 1;
                }
            }

            for (int r = 1; r < num_proc; r++)
            {
                MPI_Send(&sign, 1, MPI_INT, r, 2, MPI_COMM_WORLD);
            }
            


        } while (sign == 1);

        gettimeofday(&t_end, NULL);
        cout << "super time cost: "
            << 1000 * (t_end.tv_sec - t_start.tv_sec) +
            0.001 * (t_end.tv_usec - t_start.tv_usec) << "ms" << endl;
    }
    else
    {
        int sign;

        do
        {
            for (int i = rank; i < pasNum; i += num_proc)
            {
                MPI_Recv(&Pas[i], Num + 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            super(rank, num_proc);
            for (int i = rank; i < pasNum; i += num_proc)
            {
                MPI_Send(&Pas[i], Num + 1, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
            }

            MPI_Recv(&sign, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        } while (sign == 1);
    }
}




int main()
{
    init_A();
    init_P();
    f_ordinary1();

    init_A();
    init_P();
    f_ordinary();


    init_A();
    init_P();
    MPI_Init(NULL, NULL);

    f_mpi();

    MPI_Finalize();

}


