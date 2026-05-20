#ifndef _CONSTANT
#define _CONSTANT

#define MAX_DIMENSION 3                                          // 空間の最大次元数 
#define SIGA_DATA_RECEPTION 1                                    // SIGA解析からのデータの受け取り,  0:変位場のみ, 1:変位場と要素の重なり情報 
#define FRACTURE_ANALYSIS 2                                      // 0 : SIF計算のみ, 1 : き裂進展解析のみ 2: 両方
#define CALC_LOAD_CYCLES 0                                       // 繰り返し荷重サイクル数の計算 0:進展前のa-N曲線の傾きを使用, 1:進展前後の傾きの平均を使用
#define USE_TOTAL true                                              // 最小二乗法によるコントロールポイント更新を0 : 全てのパッチで 1 : 各パッチごとに 行う 
#define DM 1                                                     // 平面応力状態:DM = 0, 平面ひずみ状態:DM = 1
#define NG 4                                                     // Gauss-Legendreの積分点数
#define MODE_EX 1                                                // 一部で積分点数を増やして計算を行う 0, 行わない 1
#define NG_EXTEND 10                                             // 増加時のGauss-Legendreの積分点数
#define ERROR -999                                               // エラー時の返り値
#define MERGE_ERROR 1.0e-13                                      // 重なり判定の許容誤差
#define MERGE_ERROR_PARA_COORD 1.0e-10                           // パラメータ座標計算のマージの閾値
#define M_PI 3.14159265358979323846                              // 円周率
constexpr int NG3 = (NG * NG * NG);                              // NGのDIMENSION乗の計算
constexpr int NG_EXTEND3 = (NG_EXTEND * NG_EXTEND * NG_EXTEND);  // NG_EXTENDEDDのDIMENSION乗の計算
#define MAX_POW_NG NG3                                           // NGのDIMENSION乗の最大値の計算
#define MAX_POW_NG_EXTEND NG_EXTEND3                             // NGのDIMENSION乗の最大値の計算

#endif