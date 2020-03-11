#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//// 开始定义常数
#define SK       4
#define MK       25
#define HSize    256 // 0000(4) ~ 3333(4)
#define MAXINTVL 500 // 最大区间 100 bps
#define MININTVL 3   // 最小区间 3 + SK bps
#define AAAA     0
#define CCCC     85
#define GGGG     170
#define TTTT     255 // 单碱基重复序列的哈希值
//// 结束定义常数


//// 定义数据结构 Kmer
//// 这个结构用来储存当前的 Kmer
//// 结构中包含一个 int 数组，用来储存转化后的 base 值
//// 为了节约运算资源，数组采用环形结构储存
//// Kmer 的插入位置标记储存在 ins (insert) 中
typedef struct kmer *Kmer;
struct kmer {
    int mer[SK]; // 储存 base 的数组 其中 A=0 C=1 G=2 T=3
    int i;        // 表示环形数组的插入位置
};
//// 用来创建 Kmer 结构
Kmer new_Kmer(void) {
    Kmer res;
    res = (Kmer)malloc(sizeof(struct kmer));
    res->i = 0;
    return res;
}
//// 用来释放 Kmer 结构的内存空间
void free_Kmer(Kmer k) {
    free(k); // 因为 Kmer 数据结构中没有其他的指针，只需要这样既可
}
//// 将新的位值 n 添加到 Kmer 的结尾处
void push_Kmer(Kmer k, int n) {
    k->mer[k->i] = n; // 在插入位置插入 n
    k->i = (k->i + 1) % SK; // 推进插入位置
}
//// 用来将 Kmer 哈希到表值
int hash_Kmer(Kmer k) {
    int i, j, res;
    j = 0;
    res = 0;
    for ( i = k->i ; i < k->i + SK ; i ++ ) {
        res += pow(4, j) * k->mer[i%SK];
        j ++;
    }
    return res;
}
//// 完成定义数据结构 Kmer



//// 定义数据结构 PairsTab
//// 先定义 pairs 记录的行，这个行有两个值和一个指向下一个记录的指针
typedef struct pairstab_record* PairsTab_record;
struct pairstab_record {
    long pos;    // 第二个 MH 的左端位置
                 // 第二个是为了方便之后 prior 推算更长长度 MH pair 的情况需要按顺序排列
    int len;     // MH 的长度
    int ind;     // INDEL 的事件长度 即两个MH左端碱基之间的距离
    int is_mono; // 是否为单核苷酸重复
    PairsTab_record next;
};
//// 用来新建 PairsTab 记录
PairsTab_record new_PairsTab_record(long p, int l, int i, int m) {
    PairsTab_record res;
    res = (PairsTab_record)malloc(sizeof(struct pairstab_record));
    res->pos = p; res->len = l; res->ind = i;; res->is_mono = m;
    res->next = NULL; // 次项指针初始值为 NULL
    return res;
}
//// 再定义 Pairs 表格数据结构
//// 这个表格需要经历线性查询
typedef struct pairstab* PairsTab;
struct pairstab {
    PairsTab_record rec;
};
//// 用来创建 PairsTab 结构
PairsTab new_PairsTab(void) {
    PairsTab res;
    res = (PairsTab)malloc(sizeof(struct pairstab));
    res->rec = NULL;
    return res;
}
//// 用来释放 PairsTab 结构的内存空间
void free_PairsTab(PairsTab pt) {
    PairsTab_record shift;
    while (pt->rec) {
        shift = pt->rec; // 保存表头的位置到 shift
        pt->rec = pt->rec->next; // 移动表头到下一个位置
        free(shift); // 释放 shift 指向的空间
    }
    free(pt); // 释放整个 PairsTab 结构
}
//// 用来向链表的头部插入记录
void append_PairsTab(PairsTab pt, PairsTab_record r) {
    r->next = pt->rec;
    pt->rec = r;
}
//// 表格中的记录中 第二个 MH 的起始位置是逆序的
//// 法扩展 PairsTab 的记录到更长的 MH
//// 这首先要定义一个函数用来合并当前记录与相邻记录
PairsTab_record merge_PairsTab_record(PairsTab_record r) {
    int ind = r->ind;  // 合并的条件之一是 indel size 相等
    int len = SK;      // 合并最初的 MH length 为 SK
    long pos = r->pos; // 最初的 pos 为当前记录的位置
    int is_mono = r->is_mono;
    PairsTab_record res;
    for ( r = r->next ; r && (r->pos >= pos - 1 || r->pos == 0) && len < ind + SK - 1 ; r = r->next ) {
        if (r->pos == pos-1 && r->ind == ind) {
            pos--;
            len++;
            is_mono = is_mono && r->is_mono; // 如果合并的相邻 pair 全都是单碱基重复的话 合并的 pair 一定也是
            r->pos = 0;
        }
    }
    res = new_PairsTab_record(pos, len, ind, is_mono);
    return res; // 如果没合并的话，返回的记录为原记录的 copy
}
void expand_PairsTab(PairsTab pt) {
    PairsTab res = new_PairsTab(); // 新建一个栈 res 用来转移结果
    PairsTab_record holder;
    PairsTab_record merged;
    while (pt->rec) { // 从 pt 的第一项开始 直到最后一项
        if (pt->rec->pos) { // 如果当前记录没有被合并  以 pos == 0 标记
            merged = merge_PairsTab_record(pt->rec);
            if ((merged->is_mono) || (merged->len > MK) || (merged->ind - merged->len < MININTVL)) {
               // 如果是单碱基重复  或者  MH序列过长        或者        间隔序列过短的话      
                free(merged);                 // 释放内存
            } else {
                append_PairsTab(res, merged); // 将合并好的记录插入 res
            }
        }
        holder = pt->rec;
        pt->rec = pt->rec->next;
        free(holder); // 从 pt 移除当前记录
    }
    pt->rec = res->rec;
}
//// 完成定义数据结构 PairsTab




//// 定义数据结构 KHash
//// KHash 的键值为 k-mer 的4进制数值
//// KHash 的表值为 k-mer 出现的位置
//// 先定义储存表值的数据结构，这个结构用一个链表来实现
typedef struct khash_value* KHash_value;
struct khash_value {
    long pos;          // k-mer 的位置
    KHash_value next;  // 链表下一项的指针
};
//// 用来新建表值元素
KHash_value new_KHash_value(long p) {
    KHash_value res;
    res = (KHash_value)malloc(sizeof(struct khash_value));
    res->pos = p;
    res->next = NULL; // next 初始值为 NULL
    return res;
}
//// 释放哈希表值链表某个节点后的全部内存
void free_KHash_value(KHash_value pos) {
    KHash_value cut = pos;
    while (pos) {
        pos = pos->next;
        free(cut);
        cut = pos;
    }
}
//// 再定义哈希表结构
typedef struct khash* KHash;
struct khash {
    KHash_value pos[HSize];
};
//// 新建哈希表 KHash，并初始化每一个表值元素链表
KHash new_KHash(void) {
    KHash res;
    res = (KHash)malloc(sizeof(struct khash));
    int i;
    for ( i = 0 ; i < HSize ; i ++ ) {
        res->pos[i] = NULL; // 哈希表值初始化为 NULL
    }
    return res;
}
//// 用来释放哈希表的内存
void free_KHash(KHash h) {
    int i;
    for ( i = 0 ; i < HSize ; i++ ) {
        free_KHash_value(h->pos[i]);
    }
    free(h);
    //// 注意此函数只清空了指针指向的内存 并不负责将指针重新赋值与 NULL
}
//// 检查哈希表 h 的 v 位置，移除超过范围的记录
//// 并将生成的组合插入表 pt，同时在 h 的 v 位置添加新记录 p
void add_KHash_PairsTab(KHash h, int v, long p, PairsTab pt) {
    //// 将链表的顺序颠倒   新 pos 值插入在表头
    //// 这样向后扫描一遍即可 检查可否构成组合
    //// 从第一个不可构成组合的 pos 开始全部删除
    KHash_value new, shift;
    PairsTab_record holder;
    int is_mono;

    new = new_KHash_value(p); // 在表头插入新节点
    new->next = h->pos[v];
    h->pos[v] = new; 

    shift = h->pos[v]; // 拷贝链表头到 shift  此时表头为当前 pos
    while (shift->next) { // 直到 shift 走到表尾
        if ( p - shift->next->pos > MAXINTVL ) {         // 如果下一项已经超出最大检测范围
            free_KHash_value(shift->next);               // 那么从这个 pos 之后一定都超出范围了 将他们全部释放
            shift->next = NULL;                          // 并重新标记表尾的 NULL 指针
        } else {
            if ( p - shift->next->pos >= MININTVL+SK ) { // 那么在没有超过范围且大于最小界限的情况下
                if (v == AAAA || v == CCCC || v == GGGG || v == TTTT) {
                    is_mono = 1;
                } else {
                    is_mono = 0;
                }
                append_PairsTab(pt, new_PairsTab_record(p, SK, p - shift->next->pos, is_mono));
            }                                            // 向 PairsTab 添加新的记录
            shift = shift->next;                         // 向前移动 shift
        }
    }
}
//// 完成定义数据结构 KHash


//// 用来 debug
// void print_Kmer(Kmer k, KHash h) {
//     int i;
//     for ( i = k->i ; i < k->i + SK ; i++ ) {
//         switch(k->mer[i%SK]) {
//             case 0 : putchar('A') ;break;
//             case 1 : putchar('C'); break;
//             case 2 : putchar('G'); break;
//             case 3 : putchar('T'); break;
//         }
//     }
//     printf(" HASH:%5d | ", hash_Kmer(k));
//     KHash_value v = h->pos[hash_Kmer(k)];
//     while (v) {
//         printf("%5ld", v->pos);
//         v = v->next;
//     }
//     putchar('\n');
// }
//// 用来输出结果
void print_PairsTab(PairsTab pt, FILE * f , char sname[255] ) {
    PairsTab_record r;
    for ( r = pt->rec ; r ; r = r->next ) {
        fprintf(f, "%ld\t%ld\t%d\t%s\n", r->pos - r->ind, r->pos, r->len , sname );
        // printf("%10ld%10d%10d\n", r->pos, r->ind, r->len); // debug
    }
}



//// 主函数
int main(int argc, char * argv[]) {
    int c, i, j, e;        // 寄存单个字符的 c 循环变量 i j e
    int ns = 0;            // 当前序列的个数
    long pos;              // 当前序列的位置
    long skip;
    char sname[255];       // 序列的名字 Sequence NAME
    char fnout[255];          // 输出文件名 File Name OUT
    FILE * fout;           // 输出文件的指针
    Kmer k = new_Kmer();
    KHash h;
    PairsTab pt;           // 实现算法的主要数据结构

    if (argc != 2) {
        printf("Usage: find_mh <output prefix>\n");
        // 提示用法     find_mh 输出文件的前缀
    } else {
        while (1) {
            c = getchar(); // 从标准输入抓取一个字符
            switch(c) {
                case EOF :                          // 输入流结束 或者
                case '>' :                          // 检测到标记序列开始的符号 >
                    if (ns > 0) {                   // 如果 序列数已经大于 1 那么
                        // printf("%s\n", sname); // debug
                        expand_PairsTab(pt);        // 用 Priori 法扩张 PairsTab
                        print_PairsTab(pt, fout , sname);   // 输出 PairsTab
                        fclose(fout);
                        free_KHash(h);              // 释放 k-mer 哈希表
                        free_PairsTab(pt);          // 释放 PairsTab
                    }
                    if (c != EOF) {
                        scanf("%[^\n]s", sname);    // 读取新的序列名称
                        ns ++;                      // 记录过的序列数 + 1
                        i = 0;
                        while (sname[i]) {
                            if (sname[i] == ' ' || sname[i] == '\t') {
                                sname[i] = '\0';
                                break;
                            }
                            i++;
                        }                           // 在第一个分隔符的位置截断序列名称
                        sprintf(fnout, "%s.%05d.%s.out", argv[1], ns, sname);
                        fout = fopen(fnout, "w");   // 打开输出文件
                        pos = 1-SK;                 // 初始化 pos 为
                                                    // 这样当读完第一个 kmer 时 pos 为 1
                        skip = 1-SK;
                        h = new_KHash();
                        pt = new_PairsTab();
                    }
                    break;
                case 'A' :                   // 除此之外如果检测到 ACGT
                    push_Kmer(k, 0); goto _base_; // 先将碱基对应的数字插入到 Kmer 中
                case 'C' :
                    push_Kmer(k, 1); goto _base_;
                case 'G' :
                    push_Kmer(k, 2); goto _base_;
                case 'T' :
                    push_Kmer(k, 3); goto _base_;
                _base_ :
                    pos++;                   // 并将 pos 增加一位 成为当前 kmer 左端的 pos
                    skip++;
                    if (skip > 0) {          // 当 pos > 0 时 说明已经读完了整个 SK 长的 kmer
                        // print_Kmer(k, h);
                        add_KHash_PairsTab(h, hash_Kmer(k), pos, pt);
                        // print_Kmer(k, h);
                                             // 此时用计算出 kmer 的 hash 在哈希表 h 中添加 pos
                                             // 同时将满足条件的组合添加到 PairsTab pt
                    }
                    break;
                case 'N' : pos++; skip = 1-SK; break;
                default : break;
            }
            if (c == EOF) {
                break;
            }
        }
        free_Kmer(k);
    }
   
    return 0;
}
