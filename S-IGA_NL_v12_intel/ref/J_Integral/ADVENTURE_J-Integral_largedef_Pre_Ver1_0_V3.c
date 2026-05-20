/**************************************************************************************************
2016/12 荒井修正版   元のコード”ADVENTURE_J-Integral_largedef_ver5_1.c”
1)要素タイプH8,H20,T4,T10に対応。※※※T10とH8しか検証してません！！！！！！！！T4は局所最小二乗法の逆行列演算でdet=0になる場合があり、そこで死ぬ。（未修正）※※※
2)過去の拡張の残骸と思われる不要な処理を片っ端から削除。
3)カテゴリ別に関数の整理
4)内部処理が古のInHouseコード用の要素コネクティビティ、積分点の位置で行われていたが
.msh形式のコネクティビティ＋ADVENTURE出力結果の積分点位置で統一した。

argv[1]:.msh
argv[2]:.fgr
argv[3]:MARC_FEMV02.1で作成したFEMの解析結果をまとめたファイル(_JI_InputFEMResultData.dat)
argv[4]:0 or 1 (0.配置のとり方を初期配置 or 1.現在配置 にする。初期配置は節点座標そのまま。現在配置は節点座標に変位量を足しこむ)←旧越間コード、formulation checkを統合しました。
argv[5]:1 or 2 (1.微小変形問題 or 2.大変形問題)
argv[6]:1 or 2 or 3 (ガウス点の情報を各節点にふりわける方法の選択。1.局所最小二乗（線形近似式使用), 2.局所最小二乗法(二次近似式使用), 3.通常の平均計算使用)
argv[7]:Crack Flag(各節点がき裂面の上下どちらにあるのかを判定するフラグ)
argv[8]:J_DI_input_data	(temp_J_int.dat)
argv[9]:max RR(最大の積分半径)

*************************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define ERROR               -999
#define MAX_NO_NODES_ON_ELEMENT   10 //1要素あたりの最大節点数(H8※H20が死んでるのでH8で仮置き（省メモリ化のため)
#define NO_NODES_ON_ELEMENT_T10   10
#define NO_NODES_ON_ELEMENT_T4   4
#define NO_NODES_ON_ELEMENT_H8   8
#define NO_NODES_ON_ELEMENT_H20   20
#define DIMENSION           3
#define KIEL_SIZE_T10      10*DIMENSION
#define KIEL_SIZE_T4		4*DIMENSION
#define KIEL_SIZE_H8      8*DIMENSION
#define KIEL_SIZE_H20	20*DIMENSION

#define D_MATRIX_SIZE       6

#define N_STRAIN            6
#define N_STRESS            6
#define MAX_N_ELEMENT       3000000
#define MAX_N_NODE          3000000
#define MAX_LIST_ELEMENT       200000
#define MAX_LIST_NODE          2500000
#define MAX_LIST_CRACK_NODE    10000

#define MAX_N_SP 10000  //最大サンプリングポイント数
//#define N_COMPONENT 9 //サンプリングポイントが持つ物理量の数
#define MAX_N_PATCH 1000000 //最大パッチ数
#define N_P 4 //近似多項式の基底の数
#define N_P_2 10 //近似多項式の基底の数（二次基底）
#define DISPLAY 10000 //読み込み・書き出しの際にn要素(もしくはn節点)毎に処理状況をterminalに表示する。必要に応じて変更するように。
#define MAX_CRACKS 20

/*read fgr*/
#define MAX_GROUPS 40
#define MAX_FACES 100000
#define MAX_FACE_NODES 8
#define MAX_INTEGRATION_POINT 27

static int Number_of_Nodes_on_Elements;
static int N_Face_Nodes;
static int ElementShape=0;
static int ElementOrder;
static int N_IntegrationPoint;
static int N_Vertex_Nodes;

/*read msh*/
static int NElements,NNodes;
static int ElementNodeId_s[MAX_N_ELEMENT][MAX_NO_NODES_ON_ELEMENT];
static double Node_Coordinate[MAX_N_NODE][DIMENSION];

static int DispMode;//X+D or X Mode Select.

/*read InputData From FEM*/
static double Displacement[MAX_N_NODE][DIMENSION];
static double Strain[MAX_N_ELEMENT][MAX_INTEGRATION_POINT][N_STRAIN];
//越間変更点//
static double Stress[MAX_N_ELEMENT][MAX_INTEGRATION_POINT][DIMENSION][DIMENSION];
static double Pai[MAX_N_ELEMENT][MAX_INTEGRATION_POINT][DIMENSION][DIMENSION];
static double Def_grad[MAX_N_ELEMENT][MAX_INTEGRATION_POINT][DIMENSION][DIMENSION]={0};
static double Def_grad_In[MAX_N_ELEMENT][MAX_INTEGRATION_POINT][DIMENSION][DIMENSION];
static double I[MAX_N_ELEMENT][MAX_INTEGRATION_POINT][DIMENSION][DIMENSION];
//冨部変更点//
static int NElements_DI=0; //積分領域の要素総数
static int List_Element_DI[MAX_LIST_ELEMENT];
static int NNodes_DI=0; //積分領域の節点総数
static int List_Node_DI[MAX_LIST_NODE];
static int CrackNodes_DI=0; //き裂前縁の総節点数
static int List_CrackNode_DI[MAX_LIST_CRACK_NODE];
static int NCrack_Front=0;	//Number of crack front
static int Crack_N_Segments[MAX_CRACKS]={0};	// Number of Segments on Each Crack
static double Max_RR;
static double av_d;
static double RRR;
//static double Stress[MAX_N_ELEMENT][MAX_INTEGRATION_POINT][N_STRESS];
static double W[MAX_N_ELEMENT][MAX_INTEGRATION_POINT];
static double JJ[MAX_N_ELEMENT][MAX_INTEGRATION_POINT];
static double JW[MAX_N_ELEMENT][MAX_INTEGRATION_POINT];
static double Gxi_coord[MAX_N_ELEMENT][MAX_INTEGRATION_POINT][DIMENSION];
static double Coord_of_SP[MAX_N_SP][DIMENSION];
static double Val_of_SP[MAX_N_SP][9]; //各サンプリングポイントが持つ物理量
static double Coord_of_PatchCenter[MAX_N_PATCH][DIMENSION]; //パッチ中心点の座標
static double JW_on_Node[MAX_N_NODE][1];
static double Pai_on_Node[MAX_N_NODE][9];
static double Disp_grad_1_on_Node[MAX_N_NODE][3];
static double Disp_grad_2_on_Node[MAX_N_NODE][3];
static double Disp_grad_3_on_Node[MAX_N_NODE][3];
static double W_1_on_Node[MAX_N_NODE][3];
static double W_2_on_Node[MAX_N_NODE][3];
static double W_3_on_Node[MAX_N_NODE][3];

static int Crack_Front_Node[10000];
static int Total_Crack_Front_Node;
static int SP_Element_ID[10000];

//J_Integral
static double Disp_grad_1[MAX_N_ELEMENT][MAX_INTEGRATION_POINT][DIMENSION]={0};
static double Disp_grad_2[MAX_N_ELEMENT][MAX_INTEGRATION_POINT][DIMENSION]={0};
static double Disp_grad_3[MAX_N_ELEMENT][MAX_INTEGRATION_POINT][DIMENSION]={0};
static double W_1[MAX_N_ELEMENT][MAX_INTEGRATION_POINT][DIMENSION]={0};
static double W_2[MAX_N_ELEMENT][MAX_INTEGRATION_POINT][DIMENSION]={0};
static double W_3[MAX_N_ELEMENT][MAX_INTEGRATION_POINT][DIMENSION]={0};
static double W_K_D1[MAX_N_ELEMENT][MAX_INTEGRATION_POINT][DIMENSION];
static double W_K_D2[MAX_N_ELEMENT][MAX_INTEGRATION_POINT][DIMENSION];
static double W_K_D3[MAX_N_ELEMENT][MAX_INTEGRATION_POINT][DIMENSION];

static int Crack_flag[MAX_N_NODE]={0};
static int Gauss_flag[MAX_N_ELEMENT]={0};
static int crack_node[MAX_N_NODE], Ncrack_node;
//crack front local coordinate vectors at node
static double crack_node_vec[MAX_N_NODE][3][3];

int deform, method;

#define MAX_SEGS_CRACK 5000
#define MAX_CRACKS 20

//Number of crack front
static int CrackFrontSegments_For_Jint[MAX_CRACKS][MAX_SEGS_CRACK][3]={0};
// Crack front segments (nodes) for J integral computation  NCrack
static int Crack_Group[MAX_CRACKS]={0};

//Protorype//
int Make_Def_grad_In();
int Output_Def_grad_Data();
int Output_Def_grad_In_Data();
int Output_Disp_grad_Data();
int Make_Gradient_Matrix_Tetra_10(double B_1[DIMENSION][KIEL_SIZE_T10], double B_2[DIMENSION][KIEL_SIZE_T10], double B_3[DIMENSION][KIEL_SIZE_T10],
				 double Local_coord[DIMENSION], double X[NO_NODES_ON_ELEMENT_T10][DIMENSION], double *J );
int Make_Gradient_Matrix_Tetra_4(double B_1[DIMENSION][KIEL_SIZE_T4], double B_2[DIMENSION][KIEL_SIZE_T4], double B_3[DIMENSION][KIEL_SIZE_T4],
				 double Local_coord[DIMENSION], double X[NO_NODES_ON_ELEMENT_T4][DIMENSION], double *J );
int Make_Gradient_Matrix_Hexa_8(double B_1[DIMENSION][KIEL_SIZE_H8], double B_2[DIMENSION][KIEL_SIZE_H8], double B_3[DIMENSION][KIEL_SIZE_H8],
				 double Local_coord[DIMENSION], double X[NO_NODES_ON_ELEMENT_H8][DIMENSION], double *J );
int Make_Gradient_Matrix_Hexa_20(double B_1[DIMENSION][KIEL_SIZE_H20], double B_2[DIMENSION][KIEL_SIZE_H20], double B_3[DIMENSION][KIEL_SIZE_H20],
				 double Local_coord[DIMENSION], double X[NO_NODES_ON_ELEMENT_H20][DIMENSION], double *J );
int Get_SP_in_circle_distance(int data, int N_COMPONENT, int ip, int N_SP, int SP, double box);
int output_node_data(int data, int N_COMPONENT);

/***********************************************************************
						ファイル読み込み
***********************************************************************/
//ADVENTURE形式の.fgrファイル読み込み
static void ReadFgr (const char *fileName)
{
	FILE *fp;
	int iGroup;
	int iFace;
	int iNode,i,id;
	double Ga,Gb;

	fp = fopen (fileName, "r");
	if (fp == NULL) {
			fprintf (stderr, " file %s not found \n", fileName);
			exit (1);
	}

	fscanf (fp, "%d", &Number_of_Nodes_on_Elements);
	switch(Number_of_Nodes_on_Elements){
		case 4:{
			printf("Elemet Type=T4\n");
			N_Face_Nodes=3;
			ElementShape=3;
			ElementOrder=1;
			N_IntegrationPoint=1;
			N_Vertex_Nodes=4;
			break;
		}
		case 10:{
			printf("Elemet Type=T10\n");
			N_Face_Nodes=6;
			ElementShape=3;
			ElementOrder=2;
			N_IntegrationPoint=4;
			N_Vertex_Nodes=4;
			break;
		}
		case 8:{
			printf("Elemet Type=H8\n");
			N_Face_Nodes=4;
			ElementShape=4;
			ElementOrder=1;
			N_IntegrationPoint=8;
			N_Vertex_Nodes=8;
			break;
		}
		case 20:{
			printf("Elemet Type=H20\n");
			N_Face_Nodes=8;
			ElementShape=4;
			ElementOrder=2;
			N_IntegrationPoint=27;
			N_Vertex_Nodes=8;
			break;
		}
	}
	fclose (fp);
}

//ADVENTURE形式の.mshファイル読み込み
void read_file_msh(const char *fileName)
{
	FILE *fp;
	int iElement;
	int iNode;
	int i;

	fp = fopen (fileName, "r");
	if (fp == NULL) {
			fprintf (stderr, " file %s not found ¥n", fileName);
			exit (1);
	}

	fscanf (fp, "%d", &NElements);

	switch(Number_of_Nodes_on_Elements){
	  case 4:{
		for (iElement = 0; iElement < NElements; iElement++)
		{
		  fscanf (fp, "%d %d %d %d", &ElementNodeId_s[iElement][0],&ElementNodeId_s[iElement][1],&ElementNodeId_s[iElement][2],&ElementNodeId_s[iElement][3]);
		}
		break;
	  }
	  case 10:{
		for (iElement = 0; iElement < NElements; iElement++)
		{
		  fscanf (fp, "%d %d %d %d %d %d %d %d %d %d", &ElementNodeId_s[iElement][0],&ElementNodeId_s[iElement][1],
			&ElementNodeId_s[iElement][2],&ElementNodeId_s[iElement][3],&ElementNodeId_s[iElement][4],&ElementNodeId_s[iElement][5],
			&ElementNodeId_s[iElement][6],&ElementNodeId_s[iElement][7],&ElementNodeId_s[iElement][8],&ElementNodeId_s[iElement][9]);
		}
		break;
	  }
	  case 8:{
		for (iElement = 0; iElement < NElements; iElement++){
		  for (iNode = 0; iNode < Number_of_Nodes_on_Elements; iNode++){
			fscanf (fp, "%d", &ElementNodeId_s[iElement][iNode]);
		  }
		}
		break;
	  }
	  case 20:{
		for (iElement = 0; iElement < NElements; iElement++){
		  for (iNode = 0; iNode < Number_of_Nodes_on_Elements; iNode++){
			fscanf (fp, "%d", &ElementNodeId_s[iElement][iNode]);
		  }
		}
		break;
	  }
  }

	fscanf (fp, "%d", &NNodes);
	for (iNode = 0; iNode < NNodes; iNode++) {
		fscanf (fp, "%le %le %le",
			&Node_Coordinate[iNode][0],&Node_Coordinate[iNode][1],&Node_Coordinate[iNode][2]);
	}
	fclose(fp);
}

//MARC_FEMV02.1の出力形式で書きだした解析結果を読み込む。
void Get_InputData_from_FEM(const char *fileName){
	
	int i,j,k,l;char s[512];
	FILE *fp;

	if ((fp = fopen(fileName, "r")) == NULL)    printf("file open error!!¥n");
	//変位
	for(j = 0; j < NNodes; j++ ){
		for(i = 0; i < DIMENSION; i++ ){
			fscanf(fp, "%le ",&Displacement[j][i]);
		}
	}
	//ひずみ
	for( i = 0; i < NElements; i++ ){
		for( k = 0; k < N_IntegrationPoint; k++){
			for( j = 0; j < N_STRAIN; j++ ){
				fscanf(fp, "%le ",&Strain[i][k][j]);
			}
		}
	}
	//応力
	for( i = 0; i < NElements; i++ ){
		for( k = 0; k < N_IntegrationPoint; k++){
		  //応力を6成分から9成分に//
		  fscanf(fp, "%le ",&Stress[i][k][0][0]);
		  fscanf(fp, "%le ",&Stress[i][k][1][1]);
		  fscanf(fp, "%le ",&Stress[i][k][2][2]);
		  fscanf(fp, "%le ",&Stress[i][k][0][1]);
		  fscanf(fp, "%le ",&Stress[i][k][1][2]);
		  fscanf(fp, "%le ",&Stress[i][k][2][0]);
		  Stress[i][k][1][0]=Stress[i][k][0][1];
		  Stress[i][k][2][1]=Stress[i][k][1][2];
		  Stress[i][k][0][2]=Stress[i][k][2][0];
		}
	}
	//ひずみエネルギー密度
	for( i = 0; i < NElements; i++ ){
		for( k = 0; k < N_IntegrationPoint; k++ ){
			fscanf(fp, "%le", &W[i][k]);
		}
	}
	  fclose(fp);
}

//き裂前縁節点番号の読み込み
void reading_J_DI_input_data(const char *fileName)
{
  	int iCrack, iiSeg, iElement;
		int CrackNode,icp,icp_checker,ie_count;
  	int ii,jj,kk,iOrder;
  	int temp_order_of_E;
  	FILE *fp_crack;
  	fp_crack = fopen(fileName,"r");
  	int EOrder=ElementOrder+1;
		double local_vec;

	fscanf(fp_crack," %d",&temp_order_of_E);
  	fscanf(fp_crack," %d",&NCrack_Front);
		ie_count = 0;

  	for(iCrack=0; iCrack < NCrack_Front; iCrack++)
    	{
      	fscanf(fp_crack," %d", &Crack_N_Segments[iCrack]);
      	for(iiSeg=0; iiSeg < Crack_N_Segments[iCrack]; iiSeg++){
					ie_count ++;
			   for(iOrder=0;iOrder<EOrder;iOrder++){
	  			fscanf(fp_crack,"%d",&CrackNode);
					icp_checker=0;
					for(icp=0;icp<CrackNodes_DI+1;icp++){
						if(CrackNode==0) break;
						if(List_CrackNode_DI[icp]==CrackNode){
							icp_checker ++;
							break;
						}
					}
					if(icp_checker==0){
						List_CrackNode_DI[CrackNodes_DI]=CrackNode;
						CrackNodes_DI ++;
					}
	  		}

	  		for(ii=0; ii<EOrder; ii++)
	    		{
	      		for(jj=0; jj<3; jj++){
					for(kk=0; kk<3; kk++){
						fscanf(fp_crack,"  %lf",  &local_vec);
					}
	     	 	}
	    		}
		}
    	}
		
		printf("CrackFront = %d\n",CrackNodes_DI);
		printf("CrackFrontNodes\n");
		for(icp=0;icp<CrackNodes_DI;icp++){
			printf("%d ",List_CrackNode_DI[icp]);
		}
		printf("\n--------------------------------------------------\n");

  	double xc,yc,zc;
  	double xxc,yyc,zzc;
  	double d_x,d_y,d_z,D;
  	for(icp=0;icp<CrackNodes_DI-2;icp++){
    		xc = Node_Coordinate[List_CrackNode_DI[icp]][0];
    		yc = Node_Coordinate[List_CrackNode_DI[icp]][1];                     //き裂先端節点座標
    		zc = Node_Coordinate[List_CrackNode_DI[icp]][2];
    		xxc = Node_Coordinate[List_CrackNode_DI[icp+2]][0];
    		yyc = Node_Coordinate[List_CrackNode_DI[icp+2]][1];
    		zzc = Node_Coordinate[List_CrackNode_DI[icp+2]][2];
    		d_x=xxc-xc;
    		d_y=yyc-yc;
    		d_z=zzc-zc;
    		D+=sqrt(d_x*d_x + d_y*d_y + d_z*d_z);
  	}
  	av_d=D/(CrackNodes_DI-3);
  	printf("crack front average %10.10e\n",av_d);
		RRR = av_d*(Max_RR+5);
		printf("(RR+5)*av_d = %10.10e\n",RRR);

  	fclose(fp_crack);
}

/*make_nd_crack_side_V7_separate0to1and2_V2.xで作成した節点のき裂の上下判定したフラグを読み込む*/
void read_Crack_flag(const char *fileName)
{
	int iNode,nnode;
	FILE *fp_face;

	fp_face = fopen(fileName,"r");

	if (fp_face == NULL){
			fprintf (stderr, "file %s not found \n", fileName);
			exit(1);
	}

	fscanf(fp_face,"%d",&nnode);
	assert(NNodes==nnode);

	for(iNode = 0; iNode < NNodes; iNode++){
			fscanf(fp_face,"%d %d",&nnode,&Crack_flag[iNode]);
			//printf("nnode=%d check_node_flag=%d\n",nnode, Crack_flag[iNode]);
	}
	fclose(fp_face);
}
/*ガウス点（というか要素）にき裂の上下フラグを与える.read_Crack_flagとセット。*/
void make_Gauss_flag()
{

	int i, iElement, j, k;
	int count0, count1, count2;


	for(iElement=0; iElement < NElements; iElement++){
			count1=0;
			count2=0;

		for(i=0; i<Number_of_Nodes_on_Elements; i++){
			if(Crack_flag[ElementNodeId_s[iElement][i]]==3){
				printf ("ERROR!! node %d is flag 3\n", ElementNodeId_s[iElement][i]);
				exit(1);
			}
			if(Crack_flag[ElementNodeId_s[iElement][i]]==1){
				count1++;
			}
			else if(Crack_flag[ElementNodeId_s[iElement][i]]==2){
				count2++;
			}
		}

		if(count1>0 && count2>0){
			printf ("ERROR iElement=%d has Nodes U&L CrackFace (count1=%d count2=%d)\n", iElement, count1, count2);
			exit(1);
		}
		else if(count1>0){
				Gauss_flag[iElement]=1;
			//printf("%d %d\n", iElement, Gauss_flag[iElement]);
		}
		else if(count2>0){
				Gauss_flag[iElement]=2;
			//printf("%d %d\n", iElement, Gauss_flag[iElement]);
		}
		else{
				Gauss_flag[iElement]=3;
		}
	}
}
/***********************************************************************
					        形状関数
***********************************************************************/
/*形状関数設定T10*/
double N_Tetra_10(int I_No, double Local_coord[DIMENSION] ){
  double cash,temp;
	switch(I_No){
		case 0:	cash = Local_coord[0]*(2.0*Local_coord[0]-1.0);//2r(r-1/2)
				break;
		case 1:	cash = Local_coord[1]*(2.0*Local_coord[1]-1.0);//2s(s-1/2)
				break;
		case 2:	cash = (1.0-Local_coord[0]-Local_coord[1]-Local_coord[2])*(1.0-2.0*(Local_coord[0]+Local_coord[1]+Local_coord[2]));//2(1-r-s-t)(1/2-r-s-t)
				break;
		case 3:	cash = Local_coord[2]*(2.0*Local_coord[2]-1.0);//2t(t-1/2)
				break;
		case 4:	cash = 4.0*Local_coord[0]*Local_coord[1];//4rs
				break;
		case 5:	cash = 4.0*Local_coord[0]*(1-Local_coord[0]-Local_coord[1]-Local_coord[2]);//4r(1-r-s-t)
				break;
		case 6:	cash = 4.0*Local_coord[0]*Local_coord[2];//4rt		
				break;
		case 7:	cash = 4.0*Local_coord[1]*(1-Local_coord[0]-Local_coord[1]-Local_coord[2]);//4s(1-r-s-t)
				break;
		case 8:	cash = 4.0*Local_coord[2]*(1-Local_coord[0]-Local_coord[1]-Local_coord[2]);//4t(1-r-s-t)
				break;
		case 9:	cash = 4.0*Local_coord[1]*Local_coord[2];//4st_
				break;
		default:	cash = -999;
	}
	return cash;
}
/*形状関数の偏微分を求める。T10用I_No:msh,要素コネクティビティの順番,Coord[DIMENSION]:任意の座標値、xez:微分方向0~2*/
double dN_Tetra_10(int I_No, double Local_coord[DIMENSION], int xez){
	double cash;
	switch(I_No){
		case 0:
			if( xez == 0 )      	cash = 4.0*Local_coord[0]-1.0; //4r-1
			else               	cash = 0.0;
			break;
		case 1:
			if( xez == 1 )      	cash = 4.0*Local_coord[1]-1.0; //4s-1
			else                	cash = 0.0;
			break;
		case 2:
							  	cash = -3.0+4.0*(Local_coord[0]+Local_coord[1]+Local_coord[2]); //-3+4(r+s+t)
			break;
		case 3:
			if( xez == 2 )      	cash = 4.0*Local_coord[2]-1.0; //4t-1
			else                	cash = 0.0;
			break;
		case 4:
			if(xez == 0)			cash = 4.0*Local_coord[1]; //4s
			else if(xez == 1) 	cash = 4.0*Local_coord[0]; //4r
			else                	cash = 0.0;
			break;
		case 5:
			if( xez == 0 )      	cash = 4.0*(1.0-2.0*Local_coord[0]-Local_coord[1]-Local_coord[2]); //4*(1-2r-s-t)
			else                	cash = -4.0*Local_coord[0]; //-4r
			break;
		case 6:
			if( xez == 0 )      	cash = 4.0*Local_coord[2];//4t.
			else if( xez == 1 ) 	cash = 0.0;
			else                	cash = 4.0*Local_coord[0];//4r
			break;
		case 7:
			if( xez == 1 )      	cash = 4.0*(1.0-Local_coord[0]-2.0*Local_coord[1]-Local_coord[2]); //4(1-r-2s-t)
			else                	cash = -4.0*Local_coord[1]; //-4s
			break;
		case 8:
			if( xez == 2 )      	cash = 4.0*(1.0-Local_coord[0]-Local_coord[1]-2.0*Local_coord[2]); //4(1-r-s-2t)
			else                	cash = -4.0*Local_coord[2]; //-4t.
			break;
		case 9:
			if( xez == 0 )     	cash = 0.0;
			else if( xez == 1 ) 	cash = 4.0*Local_coord[2];//4t.
			else                	cash = 4.0*Local_coord[1];//4s
			break;
		default:cash = -999;
	}
	return cash;
}
/*形状関数設定T4*/
double N_Tetra_4(int I_No, double Local_coord[DIMENSION] ){
	double cash;
	switch(I_No){
		case 0: cash = Local_coord[0];   //r
			break;
		case 1: cash = Local_coord[1];  //s
			break;
		case 2: cash = 1.0-Local_coord[0]-Local_coord[1]-Local_coord[2]; //1-r-s-t
			break;
		case 3: cash = Local_coord[2]; //t
			break;
		default: cash = ERROR;
	}
	return cash;
}

/*形状関数の偏微分を求める。T4用I_No:msh,要素コネクティビティの順番,Coord[DIMENSION]:任意の座標値、xez:微分方向0~2*/
double dN_Tetra_4(int I_No, double Local_coord[DIMENSION], int xez){
	double cash;
	switch(I_No){
		case 0:
			if(xez == 0) cash = 1.0;
			else if(xez != 0) cash = 0.0;
			else cash = ERROR;
			break;
		case 1:
			if(xez == 1) cash = 1.0;
			else if(xez != 1) cash = 0.0;
			else cash = ERROR;
			break;
		case 2:
			cash = -1.0;
			break;			
		case 3:
			if(xez == 2) cash = 1.0;
			else if(xez != 2) cash = 0.0;
			else cash = ERROR;
			break;
	}
	return cash;
}

/*形状関数設定H8*/
double N_Hexa_8(int I_No, double Local_coord[DIMENSION]){
	double N;

	switch(I_No){
		case 0:
	  		N = 0.125 * ( 1.0 - Local_coord[0] ) * ( 1.0 - Local_coord[1] ) * ( 1.0 - Local_coord[2] ) ;
	  		break;
	  	case 1:
 			N = 0.125 * ( 1.0 + Local_coord[0] ) * ( 1.0 - Local_coord[1] ) * ( 1.0 - Local_coord[2] ) ;
 			break;
	  	case 2:
  			N = 0.125 * ( 1.0 + Local_coord[0] ) * ( 1.0 + Local_coord[1] ) * ( 1.0 - Local_coord[2] ) ;
  			break;
	  	case 3:
  			N = 0.125 * ( 1.0 - Local_coord[0] ) * ( 1.0 + Local_coord[1] ) * ( 1.0 - Local_coord[2] ) ;
  			break;
	  	case 4:
  			N = 0.125 * ( 1.0 - Local_coord[0] ) * ( 1.0 - Local_coord[1] ) * ( 1.0 + Local_coord[2] ) ;
  			break;
	  	case 5:
  			N = 0.125 * ( 1.0 + Local_coord[0] ) * ( 1.0 - Local_coord[1] ) * ( 1.0 + Local_coord[2] ) ;
  			break;
	  	case 6:
  			N = 0.125 * ( 1.0 + Local_coord[0] ) * ( 1.0 + Local_coord[1] ) * ( 1.0 + Local_coord[2] ) ;
  			break;
	  	case 7:
  			N = 0.125 * ( 1.0 - Local_coord[0] ) * ( 1.0 + Local_coord[1] ) * ( 1.0 + Local_coord[2] ) ;
  			break;
  		default: N = ERROR;
	}
	return N;
}

//形状関数の偏微分H8 (I_No:節点番号 xez:偏微分の分母部分0ξ1η2Γ）
double dN_Hexa_8 (int I_No, double Local_coord[DIMENSION],int xez){
	double dN;

	switch(I_No){
		case 0:
			switch(xez){
				case 0: dN = -0.125 * ( 1.0 - Local_coord[1] ) * ( 1.0 - Local_coord[2] ) ;
					break;
				case 1: dN = -0.125 * ( 1.0 - Local_coord[0] ) * ( 1.0 - Local_coord[2] ) ;
					break;
				case 2: dN = -0.125 * ( 1.0 - Local_coord[0] ) * ( 1.0 - Local_coord[1] ) ;
					break;
				default: dN=ERROR;
			}
			break;
		case 1:
			switch(xez){
				case 0: dN =  0.125 * ( 1.0 - Local_coord[1] ) * ( 1.0 - Local_coord[2] ) ;
					break;
				case 1: dN = -0.125 * ( 1.0 + Local_coord[0] ) * ( 1.0 - Local_coord[2] ) ;
					break;
				case 2: dN = -0.125 * ( 1.0 + Local_coord[0] ) * ( 1.0 - Local_coord[1] ) ;
					break;
				default: dN=ERROR;
			}
			break;
		case 2:
			switch(xez){
				case 0: dN =  0.125 * ( 1.0 + Local_coord[1] ) * ( 1.0 - Local_coord[2] ) ;
					break;
				case 1: dN =  0.125 * ( 1.0 + Local_coord[0] ) * ( 1.0 - Local_coord[2] ) ;
					break;
				case 2: dN = -0.125 * ( 1.0 + Local_coord[0] ) * ( 1.0 + Local_coord[1] ) ;
					break;
				default: dN=ERROR;
			}
			break;
		case 3:
			switch(xez){
				case 0: dN = -0.125 * ( 1.0 + Local_coord[1] ) * ( 1.0 - Local_coord[2] ) ;
					break;
				case 1: dN =  0.125 * ( 1.0 - Local_coord[0] ) * ( 1.0 - Local_coord[2] ) ;
					break;
				case 2: dN = -0.125 * ( 1.0 - Local_coord[0] ) * ( 1.0 + Local_coord[1] ) ;
					break;
				default: dN=ERROR;
			}
			break;
		case 4:
			switch(xez){
				case 0:  dN = -0.125 * ( 1.0 - Local_coord[1] ) * ( 1.0 + Local_coord[2] ) ;
					break;
				case 1:	dN = -0.125 * ( 1.0 - Local_coord[0] ) * ( 1.0 + Local_coord[2] ) ;
					break;
				case 2:  dN =  0.125 * ( 1.0 - Local_coord[0] ) * ( 1.0 - Local_coord[1] ) ;
					break;
				default: dN=ERROR;
			}
			break;
		case 5:
			switch(xez){
				case 0:  dN =  0.125 * ( 1.0 - Local_coord[1] ) * ( 1.0 + Local_coord[2] ) ;
					break;
				case 1:  dN = -0.125 * ( 1.0 + Local_coord[0] ) * ( 1.0 + Local_coord[2] ) ;
					break;
				case 2:  dN =  0.125 * ( 1.0 + Local_coord[0] ) * ( 1.0 - Local_coord[1] ) ;
					break;
				default: dN=ERROR;
			}
			break;
		case 6:
			switch(xez){
				case 0:  dN =  0.125 * ( 1.0 + Local_coord[1] ) * ( 1.0 + Local_coord[2] ) ;
					break;
				case 1:  dN =  0.125 * ( 1.0 + Local_coord[0] ) * ( 1.0 + Local_coord[2] ) ;
					break;
				case 2:  dN =  0.125 * ( 1.0 + Local_coord[0] ) * ( 1.0 + Local_coord[1] ) ;
					break;
				default: dN=ERROR;
			}
			break;
		case 7:
			switch(xez){
				case 0:  dN = -0.125 * ( 1.0 + Local_coord[1] ) * ( 1.0 + Local_coord[2] ) ;
					break;
				case 1:  dN =  0.125 * ( 1.0 - Local_coord[0] ) * ( 1.0 + Local_coord[2] ) ;
					break;
				case 2:  dN =  0.125 * ( 1.0 - Local_coord[0] ) * ( 1.0 + Local_coord[1] ) ;
					break;
				default: dN=ERROR;
			}
			break;
		default:dN=ERROR;

	}
	return dN;
}

/*形状関数設定H20*/
double N_Hexa_20(int I_No, double Local_coord[DIMENSION]){
	double N;

	switch(I_No){
		case 0: N = -0.125 * (1.0 - Local_coord[0]) * ( 1.0 - Local_coord[1] ) * ( 1.0 - Local_coord[2] ) * ( 2.0 + Local_coord[0] + Local_coord[1] + Local_coord[2] ) ;
			break;
		case 1: N = -0.125 * (1.0 + Local_coord[0]) * ( 1.0 - Local_coord[1] ) * ( 1.0 - Local_coord[2] ) * ( 2.0 - Local_coord[0] + Local_coord[1] + Local_coord[2] ) ;
			break;
		case 2: N = -0.125 * (1.0 + Local_coord[0]) * ( 1.0 + Local_coord[1] ) * ( 1.0 - Local_coord[2] ) * ( 2.0 - Local_coord[0] - Local_coord[1] + Local_coord[2] ) ;
			break;
		case 3: N = -0.125 * (1.0 - Local_coord[0]) * ( 1.0 + Local_coord[1] ) * ( 1.0 - Local_coord[2] ) * ( 2.0 + Local_coord[0] - Local_coord[1] + Local_coord[2] ) ;
			break;
		case 4: N = -0.125 * (1.0 - Local_coord[0]) * ( 1.0 - Local_coord[1] ) * ( 1.0 + Local_coord[2] ) * ( 2.0 + Local_coord[0] + Local_coord[1] - Local_coord[2] ) ;
			break;
		case 5: N = -0.125 * (1.0 + Local_coord[0]) * ( 1.0 - Local_coord[1] ) * ( 1.0 + Local_coord[2] ) * ( 2.0 - Local_coord[0] + Local_coord[1] - Local_coord[2] ) ;
			break;
		case 6: N = -0.125 * (1.0 + Local_coord[0]) * ( 1.0 + Local_coord[1] ) * ( 1.0 + Local_coord[2] ) * ( 2.0 - Local_coord[0] - Local_coord[1] - Local_coord[2] ) ;
			break;
		case 7: N = -0.125 * (1.0 - Local_coord[0]) * ( 1.0 + Local_coord[1] ) * ( 1.0 + Local_coord[2] ) * ( 2.0 + Local_coord[0] - Local_coord[1] - Local_coord[2] ) ;
			break;

		case 8 : N = 0.25 * (1.0 - Local_coord[0]*Local_coord[0]) * ( 1.0 - Local_coord[1] ) * ( 1.0 - Local_coord[2] ) ;
			break;
		case 9 : N = 0.25 * (1.0 + Local_coord[0]) * ( 1.0 - Local_coord[1]*Local_coord[1] ) * ( 1.0 - Local_coord[2] ) ;
			break;
		case 10: N = 0.25 * (1.0 - Local_coord[0]*Local_coord[0]) * ( 1.0 + Local_coord[1] ) * ( 1.0 - Local_coord[2] ) ;
			break;
		case 11: N = 0.25 * (1.0 - Local_coord[0]) * ( 1.0 - Local_coord[1]*Local_coord[1] ) * ( 1.0 - Local_coord[2] ) ;
			break;

		case 12: N = 0.25 * (1.0 - Local_coord[0]) * ( 1.0 - Local_coord[1] ) * ( 1.0 - Local_coord[2]*Local_coord[2] ) ;
			break;
		case 13: N = 0.25 * (1.0 + Local_coord[0]) * ( 1.0 - Local_coord[1] ) * ( 1.0 - Local_coord[2]*Local_coord[2] ) ;
			break;
		case 14: N = 0.25 * (1.0 + Local_coord[0]) * ( 1.0 + Local_coord[1] ) * ( 1.0 - Local_coord[2]*Local_coord[2] ) ;
			break;
		case 15: N = 0.25 * (1.0 - Local_coord[0]) * ( 1.0 + Local_coord[1] ) * ( 1.0 - Local_coord[2]*Local_coord[2] ) ;
			break;

		case 16: N = 0.25 * (1.0 - Local_coord[0]*Local_coord[0]) * ( 1.0 - Local_coord[1] ) * ( 1.0 + Local_coord[2] ) ;
			break;
		case 17: N = 0.25 * (1.0 + Local_coord[0]) * ( 1.0 - Local_coord[1]*Local_coord[1] ) * ( 1.0 + Local_coord[2] ) ;
			break;
		case 18: N = 0.25 * (1.0 - Local_coord[0]*Local_coord[0]) * ( 1.0 + Local_coord[1] ) * ( 1.0 + Local_coord[2] ) ;
			break;
		case 19: N = 0.25 * (1.0 - Local_coord[0]) * ( 1.0 - Local_coord[1]*Local_coord[1] ) * ( 1.0 + Local_coord[2] ) ;
			break;
		default:N=ERROR;
	}
	return N;
}
//形状関数の偏微分H20 (I_No:節点番号 xez:偏微分の分母部分0ξ1η2Γ）
double dN_Hexa_20(int I_No, double Local_coord[DIMENSION], int xez){
	double dN;

	switch(I_No)
	{
	  	case 0:
	  		switch(xez){
	  			case 0: dN = -0.125 * ( 1.0 - Local_coord[1] ) * ( 1.0 - Local_coord[2] ) * ( -2.0 * Local_coord[0] - Local_coord[1] - Local_coord[2] - 1.0 ) ;
	  				break;
	  			case 1: dN = -0.125 * ( 1.0 - Local_coord[0] ) * ( 1.0 - Local_coord[2] ) * ( -2.0 * Local_coord[1] - Local_coord[0] - Local_coord[2] - 1.0 ) ;
	  				break;
	  			case 2: dN = -0.125 * ( 1.0 - Local_coord[0] ) * ( 1.0 - Local_coord[1] ) * ( -2.0 * Local_coord[2] - Local_coord[0] - Local_coord[1] - 1.0 ) ;
	  				break;
	  			default:dN=ERROR;
	  		}
	  		break;
	  	case 1:
	  		switch(xez){
	  			case 0: dN =  0.125 * ( 1.0 - Local_coord[1] ) * ( 1.0 - Local_coord[2] ) * (  2.0 * Local_coord[0] - Local_coord[1] - Local_coord[2] - 1.0 ) ;
	  				break;
	  			case 1: dN = -0.125 * ( 1.0 + Local_coord[0] ) * ( 1.0 - Local_coord[2] ) * ( -2.0 * Local_coord[1] + Local_coord[0] - Local_coord[2] - 1.0 ) ;
	  				break;
	  			case 2: dN = -0.125 * ( 1.0 + Local_coord[0] ) * ( 1.0 - Local_coord[1] ) * ( -2.0 * Local_coord[2] + Local_coord[0] - Local_coord[1] - 1.0 ) ;
	  				break;
	  			default:dN=ERROR;
	  		}
	  		break;
	  	case 2:
	  		switch(xez){
	  			case 0: dN =  0.125 * ( 1.0 + Local_coord[1] ) * ( 1.0 - Local_coord[2] ) * (  2.0 * Local_coord[0] + Local_coord[1] - Local_coord[2] - 1.0 ) ;
	  				break;
	  			case 1: dN =  0.125 * ( 1.0 + Local_coord[0] ) * ( 1.0 - Local_coord[2] ) * (  2.0 * Local_coord[1] + Local_coord[0] - Local_coord[2] - 1.0 ) ;
	  				break;
	  			case 2: dN = -0.125 * ( 1.0 + Local_coord[0] ) * ( 1.0 + Local_coord[1] ) * ( -2.0 * Local_coord[2] + Local_coord[0] + Local_coord[1] - 1.0 ) ;
	  				break;
	  			default:dN=ERROR;
	  		}
	  		break;
	  	case 3:
	  		switch(xez){
	  			case 0: dN = -0.125 * ( 1.0 + Local_coord[1] ) * ( 1.0 - Local_coord[2] ) * ( -2.0 * Local_coord[0] + Local_coord[1] - Local_coord[2] - 1.0 ) ;
	  				break;
	  			case 1: dN =  0.125 * ( 1.0 - Local_coord[0] ) * ( 1.0 - Local_coord[2] ) * (  2.0 * Local_coord[1] - Local_coord[0] - Local_coord[2] - 1.0 ) ;
	  				break;
	  			case 2: dN = -0.125 * ( 1.0 - Local_coord[0] ) * ( 1.0 + Local_coord[1] ) * ( -2.0 * Local_coord[2] - Local_coord[0] + Local_coord[1] - 1.0 ) ;
	  				break;
	  			default:dN=ERROR;
	  		}
	  		break;
	  	case 4:
	  		switch(xez){
	  			case 0: dN = -0.125 * ( 1.0 - Local_coord[1] ) * ( 1.0 + Local_coord[2] ) * ( -2.0 * Local_coord[0] - Local_coord[1] + Local_coord[2] - 1.0 ) ;
	  				break;
	  			case 1: dN = -0.125 * ( 1.0 - Local_coord[0] ) * ( 1.0 + Local_coord[2] ) * ( -2.0 * Local_coord[1] - Local_coord[0] + Local_coord[2] - 1.0 ) ;
	  				break;
	  			case 2: dN =  0.125 * ( 1.0 - Local_coord[0] ) * ( 1.0 - Local_coord[1] ) * (  2.0 * Local_coord[2] - Local_coord[0] - Local_coord[1] - 1.0 ) ;
	  				break;
	  			default:dN=ERROR;
	  		}
	  		break;
	  	case 5:
	  		switch(xez){
	  			case 0: dN =  0.125 * ( 1.0 - Local_coord[1] ) * ( 1.0 + Local_coord[2] ) * (  2.0 * Local_coord[0] - Local_coord[1] + Local_coord[2] - 1.0 ) ;
	  				break;
	  			case 1: dN = -0.125 * ( 1.0 + Local_coord[0] ) * ( 1.0 + Local_coord[2] ) * ( -2.0 * Local_coord[1] + Local_coord[0] + Local_coord[2] - 1.0 ) ;
	  				break;
	  			case 2: dN =  0.125 * ( 1.0 + Local_coord[0] ) * ( 1.0 - Local_coord[1] ) * (  2.0 * Local_coord[2] + Local_coord[0] - Local_coord[1] - 1.0 ) ;
	  				break;
	  			default:dN=ERROR;
	  		}
	  		break;
	  	case 6:
	  		switch(xez){
	  			case 0: dN =  0.125 * ( 1.0 + Local_coord[1] ) * ( 1.0 + Local_coord[2] ) * (  2.0 * Local_coord[0] + Local_coord[1] + Local_coord[2] - 1.0 ) ;
	  				break;
	  			case 1: dN =  0.125 * ( 1.0 + Local_coord[0] ) * ( 1.0 + Local_coord[2] ) * (  2.0 * Local_coord[1] + Local_coord[0] + Local_coord[2] - 1.0 ) ;
	  				break;
	  			case 2: dN =  0.125 * ( 1.0 + Local_coord[0] ) * ( 1.0 + Local_coord[1] ) * (  2.0 * Local_coord[2] + Local_coord[0] + Local_coord[1] - 1.0 ) ;
	  				break;
	  			default:dN=ERROR;
	  		}
	  		break;
	  	case 7:
	  		switch(xez){
	  			case 0: dN = -0.125 * ( 1.0 + Local_coord[1] ) * ( 1.0 + Local_coord[2] ) * ( -2.0 * Local_coord[0] + Local_coord[1] + Local_coord[2] - 1.0 ) ;
	  				break;
	  			case 1: dN =  0.125 * ( 1.0 - Local_coord[0] ) * ( 1.0 + Local_coord[2] ) * (  2.0 * Local_coord[1] - Local_coord[0] + Local_coord[2] - 1.0 ) ;
	  				break;
	  			case 2: dN =  0.125 * ( 1.0 - Local_coord[0] ) * ( 1.0 + Local_coord[1] ) * (  2.0 * Local_coord[2] - Local_coord[0] + Local_coord[1] - 1.0 ) ;
	  				break;
	  			default:dN=ERROR;
	  		}
	  		break;

	  	case 8:
	  		switch(xez){
	  			case 0: dN = -0.5  * Local_coord[0] * (1.0 - Local_coord[1]) * (1.0 - Local_coord[2]) ;
	  				break;
	  			case 1: dN = -0.25 * (1.0 - Local_coord[0] * Local_coord[0]) * (1.0 - Local_coord[2]) ;
	  				break;
	  			case 2: dN = -0.25 * (1.0 - Local_coord[0] * Local_coord[0]) * (1.0 - Local_coord[1]) ;
	  				break;
	  			default:dN=ERROR;
	  		}
	  		break;
	  	case 9:
	  		switch(xez){
	  			case 0: dN =  0.25 * (1.0 - Local_coord[1] * Local_coord[1]) * (1.0 - Local_coord[2]) ;
	  				break;
	  			case 1: dN = -0.5  * Local_coord[1] * (1.0 + Local_coord[0]) * (1.0 - Local_coord[2]) ;
	  				break;
	  			case 2: dN = -0.25 * (1.0 - Local_coord[1] * Local_coord[1]) * (1.0 + Local_coord[0]) ;
	  				break;
	  			default:dN=ERROR;
	  		}
	  		break;
	  	case 10:
	  		switch(xez){
	  			case 0: dN = -0.5  * Local_coord[0] * (1.0 + Local_coord[1]) * (1.0 - Local_coord[2]) ;
	  				break;
	  			case 1: dN =  0.25 * (1.0 - Local_coord[0] * Local_coord[0]) * (1.0 - Local_coord[2]) ;
	  				break;
	  			case 2: dN = -0.25 * (1.0 - Local_coord[0] * Local_coord[0]) * (1.0 + Local_coord[1]) ;
	  				break;
	  			default:dN=ERROR;
	  		}
	  		break;
	  	case 11:
	  		switch(xez){
	  			case 0: dN = -0.25 * (1.0 - Local_coord[1] * Local_coord[1]) * (1.0 - Local_coord[2]) ;
	  				break;
	  			case 1: dN = -0.5  * Local_coord[1] * (1.0 - Local_coord[0]) * (1.0 - Local_coord[2]) ;
	  				break;
	  			case 2: dN = -0.25 * (1.0 - Local_coord[1] * Local_coord[1]) * (1.0 - Local_coord[0]) ;
	  				break;
	  			default:dN=ERROR;
	  		}
	  		break;
	  	case 12:
	  		switch(xez){
	  			case 0: dN = -0.25 * (1.0 - Local_coord[2] * Local_coord[2]) * (1.0 - Local_coord[1]) ;
	  				break;
	  			case 1: dN = -0.25 * (1.0 - Local_coord[2] * Local_coord[2]) * (1.0 - Local_coord[0]) ;
	  				break;
	  			case 2: dN = -0.5  * Local_coord[2] * (1.0 - Local_coord[0]) * (1.0 - Local_coord[1]) ;
	  				break;
	  			default:dN=ERROR;
	  		}
	  		break;
	  	case 13:
	  		switch(xez){
	  			case 0: dN =  0.25 * (1.0 - Local_coord[2] * Local_coord[2]) * (1.0 - Local_coord[1]) ;
	  				break;
	  			case 1: dN = -0.25 * (1.0 - Local_coord[2] * Local_coord[2]) * (1.0 + Local_coord[0]) ;
	  				break;
	  			case 2: dN = -0.5  * Local_coord[2] * (1.0 + Local_coord[0]) * (1.0 - Local_coord[1]) ;
	  				break;
	  			default:dN=ERROR;
	  		}
	  		break;
	  	case 14:
	  		switch(xez){
	  			case 0: dN =  0.25 * (1.0 - Local_coord[2] * Local_coord[2]) * (1.0 + Local_coord[1]) ;
	  				break;
	  			case 1: dN =  0.25 * (1.0 - Local_coord[2] * Local_coord[2]) * (1.0 + Local_coord[0]) ;
	  				break;
	  			case 2: dN = -0.5  * Local_coord[2] * (1.0 + Local_coord[0]) * (1.0 + Local_coord[1]) ;
	  				break;
	  			default:dN=ERROR;
	  		}
	  		break;
	  	case 15:
	  		switch(xez){
	  			case 0: dN = -0.25 * (1.0 - Local_coord[2] * Local_coord[2]) * (1.0 + Local_coord[1]) ;
	  				break;
	  			case 1: dN =  0.25 * (1.0 - Local_coord[2] * Local_coord[2]) * (1.0 - Local_coord[0]) ;
	  				break;
	  			case 2: dN = -0.5  * Local_coord[2] * (1.0 - Local_coord[0]) * (1.0 + Local_coord[1]) ;
	  				break;
	  			default:dN=ERROR;
	  		}
	  		break;
	  	case 16:
	  		switch(xez){
	  			case 0: dN = -0.5  * Local_coord[0] * (1.0 - Local_coord[1]) * (1.0 + Local_coord[2]) ;
	  				break;
	  			case 1: dN = -0.25 * (1.0 - Local_coord[0] * Local_coord[0]) * (1.0 + Local_coord[2]) ;
	  				break;
	  			case 2: dN =  0.25 * (1.0 - Local_coord[0] * Local_coord[0]) * (1.0 - Local_coord[1]) ;
	  				break;
	  			default:dN=ERROR;
	  		}
	  		break;
	  	case 17:
	  		switch(xez){
	  			case 0: dN =  0.25 * (1.0 - Local_coord[1] * Local_coord[1]) * (1.0 + Local_coord[2]) ;
	  				break;
	  			case 1: dN = -0.5  * Local_coord[1] * (1.0 + Local_coord[0]) * (1.0 + Local_coord[2]) ;
	  				break;
	  			case 2: dN =  0.25 * (1.0 - Local_coord[1] * Local_coord[1]) * (1.0 + Local_coord[0]) ;
	  				break;
	  			default:dN=ERROR;
	  		}
	  		break;
	  	case 18:
	  		switch(xez){
	  			case 0: dN= -0.5  * Local_coord[0] * (1.0 + Local_coord[1]) * (1.0 + Local_coord[2]) ;
	  				break;
	  			case 1: dN =  0.25 * (1.0 - Local_coord[0] * Local_coord[0]) * (1.0 + Local_coord[2]) ;
	  				break;
	  			case 2: dN =  0.25 * (1.0 - Local_coord[0] * Local_coord[0]) * (1.0 + Local_coord[1]) ;
	  				break;
	  			default:dN=ERROR;
	  		}
	  		break;
	  	case 19:
	  		switch(xez){
	  			case 0: dN = -0.25 * (1.0 - Local_coord[1] * Local_coord[1]) * (1.0 + Local_coord[2]) ;
	  				break;
	  			case 1: dN = -0.5  * Local_coord[1] * (1.0 - Local_coord[0]) * (1.0 + Local_coord[2]) ;
	  				break;
	  			case 2: dN =  0.25 * (1.0 - Local_coord[1] * Local_coord[1]) * (1.0 - Local_coord[0]) ;
	  				break;
	  			default:dN=ERROR;
	  		}
	  		break;
	  	default:dN=ERROR;
	}
	return dN;
}


/***********************************************************************
					     ヤコビMatrix
***********************************************************************/
int Jacobian_Tetra_10( double a[DIMENSION][DIMENSION], double Local_coord[DIMENSION], double X[NO_NODES_ON_ELEMENT_T10][DIMENSION] ){
	int i,j,k;
	for( i= 0; i < DIMENSION; i++ ){
		for(j = 0; j < DIMENSION; j++ ){
			a[i][j] = 0.0;
			for( k = 0; k < Number_of_Nodes_on_Elements; k++ ){
				a[i][j] += dN_Tetra_10( k, Local_coord, j)*X[k][i];
			}
			//printf("\na[%d][%d]=%5.25lf\n",i,j,a[i][j]);
		}
	}
	return 0;
}
int Jacobian_Tetra_4( double a[DIMENSION][DIMENSION], double Local_coord[DIMENSION], double X[NO_NODES_ON_ELEMENT_T4][DIMENSION] ){
	int i,j,k;
	for( i= 0; i < DIMENSION; i++ ){
		for(j = 0; j < DIMENSION; j++ ){
			a[i][j] = 0.0;
			for( k = 0; k < Number_of_Nodes_on_Elements; k++ ){
				a[i][j] += dN_Tetra_4( k, Local_coord, j)*X[k][i];
			}
			//printf("a[%d][%d]=%lf\n",i,j,a[i][j]);
		}
	}
	return 0;
}
int Jacobian_Hexa_8( double a[DIMENSION][DIMENSION], double Local_coord[DIMENSION], double X[NO_NODES_ON_ELEMENT_H8][DIMENSION] ){
	int i,j,k;
	for( i= 0; i < DIMENSION; i++ ){
		for(j = 0; j < DIMENSION; j++ ){
			a[i][j] = 0.0;
			for( k = 0; k < Number_of_Nodes_on_Elements; k++ ){
				a[i][j] += dN_Hexa_8( k, Local_coord, j)*X[k][i];
			}
			//printf("a[%d][%d]=%lf\n",i,j,a[i][j]);
			//printf("a[2][2]: %lf\n", a[2][2]);

		}
	}
	return 0;
}
int Jacobian_Hexa_20( double a[DIMENSION][DIMENSION], double Local_coord[DIMENSION], double X[NO_NODES_ON_ELEMENT_H20][DIMENSION] ){
	int i,j,k;
	for( i= 0; i < DIMENSION; i++ ){
		for(j = 0; j < DIMENSION; j++ ){
			a[i][j] = 0.0;
			for( k = 0; k < Number_of_Nodes_on_Elements; k++ ){
				a[i][j] += dN_Hexa_20( k, Local_coord, j)*X[k][i];
			}
			//printf("a[%d][%d]=%lf\n",i,j,a[i][j]);
		}
	}
	return 0;
}

/***********************************************************************
					     行列演算
***********************************************************************/
/*読み込んだ3×3の行列を逆行列にして返す。返り値がdetになってる。*/
double InverseMatrix_3D(double M[3][3] ){
  	int i, j;
  	double a[3][3];
  	double det = M[0][0]*M[1][1]*M[2][2] +M[0][1]*M[1][2]*M[2][0] +M[0][2]*M[1][0]*M[2][1]
				-M[0][0]*M[1][2]*M[2][1] -M[0][1]*M[1][0]*M[2][2] -M[0][2]*M[1][1]*M[2][0];
	//printf("Jacobian determinant det = %lf\n", det);

  	for( i= 0; i< 3; i++ ){
		for( j= 0; j< 3; j++ )  a[i][j] = M[i][j];
  	}
  	M[0][0] = (a[1][1]*a[2][2]-a[1][2]*a[2][1])/det; M[0][1] = (a[2][1]*a[0][2]-a[2][2]*a[0][1])/det; M[0][2] = (a[0][1]*a[1][2]-a[0][2]*a[1][1])/det;
  	M[1][0] = (a[1][2]*a[2][0]-a[1][0]*a[2][2])/det; M[1][1] = (a[2][2]*a[0][0]-a[2][0]*a[0][2])/det; M[1][2] = (a[0][2]*a[1][0]-a[0][0]*a[1][2])/det;
  	M[2][0] = (a[1][0]*a[2][1]-a[1][1]*a[2][0])/det; M[2][1] = (a[2][0]*a[0][1]-a[2][1]*a[0][0])/det; M[2][2] = (a[0][0]*a[1][1]-a[0][1]*a[1][0])/det;
  	return det;
}

/*読み込んだ4×4の行列を逆行列にして返す。返り値がdetになってる。*/
double InverseMatrix_4D( double M4D[4][4], int ip ){

  	int i, j;
  	long double a[4][4];
  	double det =   M4D[0][0]*M4D[1][1]*M4D[2][2]*M4D[3][3] + M4D[0][0]*M4D[1][2]*M4D[2][3]*M4D[3][1] + M4D[0][0]*M4D[1][3]*M4D[2][1]*M4D[3][2]
				+ M4D[0][1]*M4D[1][0]*M4D[2][3]*M4D[3][2] + M4D[0][1]*M4D[1][2]*M4D[2][0]*M4D[3][3] + M4D[0][1]*M4D[1][3]*M4D[2][2]*M4D[3][0]
				+ M4D[0][2]*M4D[1][0]*M4D[2][1]*M4D[3][3] + M4D[0][2]*M4D[1][1]*M4D[2][3]*M4D[3][0] + M4D[0][2]*M4D[1][3]*M4D[2][0]*M4D[3][1]
				+ M4D[0][3]*M4D[1][0]*M4D[2][2]*M4D[3][1] + M4D[0][3]*M4D[1][1]*M4D[2][0]*M4D[3][2] + M4D[0][3]*M4D[1][2]*M4D[2][1]*M4D[3][0]
				- M4D[0][0]*M4D[1][1]*M4D[2][3]*M4D[3][2] - M4D[0][0]*M4D[1][2]*M4D[2][1]*M4D[3][3] - M4D[0][0]*M4D[1][3]*M4D[2][2]*M4D[3][1]
				- M4D[0][1]*M4D[1][0]*M4D[2][2]*M4D[3][3] - M4D[0][1]*M4D[1][2]*M4D[2][3]*M4D[3][0] - M4D[0][1]*M4D[1][3]*M4D[2][0]*M4D[3][2]
				- M4D[0][2]*M4D[1][0]*M4D[2][3]*M4D[3][1] - M4D[0][2]*M4D[1][1]*M4D[2][0]*M4D[3][3] - M4D[0][2]*M4D[1][3]*M4D[2][1]*M4D[3][0]
				- M4D[0][3]*M4D[1][0]*M4D[2][1]*M4D[3][2] - M4D[0][3]*M4D[1][1]*M4D[2][2]*M4D[3][0] - M4D[0][3]*M4D[1][2]*M4D[2][0]*M4D[3][1];
			//if(det<0.000000001)	printf("det=%le\n",det );
  	for( i= 0 ; i< 4 ; i++ ){
		for( j= 0 ; j< 4 ; j++ ){
	  		a[i][j] = M4D[i][j];
		}
  	}
  	M4D[0][0] = ( a[1][1]*a[2][2]*a[3][3] + a[1][2]*a[2][3]*a[3][1] + a[1][3]*a[2][1]*a[3][2] - a[1][1]*a[2][3]*a[3][2] - a[1][2]*a[2][1]*a[3][3] - a[1][3]*a[2][2]*a[3][1] ) / det;
  	M4D[0][1] = ( a[0][1]*a[2][3]*a[3][2] + a[0][2]*a[2][1]*a[3][3] + a[0][3]*a[2][2]*a[3][1] - a[0][1]*a[2][2]*a[3][3] - a[0][2]*a[2][3]*a[3][1] - a[0][3]*a[2][1]*a[3][2] ) / det;
  	M4D[0][2] = ( a[0][1]*a[1][2]*a[3][3] + a[0][2]*a[1][3]*a[3][1] + a[0][3]*a[1][1]*a[3][2] - a[0][1]*a[1][3]*a[3][2] - a[0][2]*a[1][1]*a[3][3] - a[0][3]*a[1][2]*a[3][1] ) / det;
  	M4D[0][3] = ( a[0][1]*a[1][3]*a[2][2] + a[0][2]*a[1][1]*a[2][3] + a[0][3]*a[1][2]*a[2][1] - a[0][1]*a[1][2]*a[2][3] - a[0][2]*a[1][3]*a[2][1] - a[0][3]*a[1][1]*a[2][2] ) / det;
  	M4D[1][0] = ( a[1][0]*a[2][3]*a[3][2] + a[1][2]*a[2][0]*a[3][3] + a[1][3]*a[2][2]*a[3][0] - a[1][0]*a[2][2]*a[3][3] - a[1][2]*a[2][3]*a[3][0] - a[1][3]*a[2][0]*a[3][2] ) / det;
  	M4D[1][1] = ( a[0][0]*a[2][2]*a[3][3] + a[0][2]*a[2][3]*a[3][0] + a[0][3]*a[2][0]*a[3][2] - a[0][0]*a[2][3]*a[3][2] - a[0][2]*a[2][0]*a[3][3] - a[0][3]*a[2][2]*a[3][0] ) / det;
  	M4D[1][2] = ( a[0][0]*a[1][3]*a[3][2] + a[0][2]*a[1][0]*a[3][3] + a[0][3]*a[1][2]*a[3][0] - a[0][0]*a[1][2]*a[3][3] - a[0][2]*a[1][3]*a[3][0] - a[0][3]*a[1][0]*a[3][2] ) / det;
  	M4D[1][3] = ( a[0][0]*a[1][2]*a[2][3] + a[0][2]*a[1][3]*a[2][0] + a[0][3]*a[1][0]*a[2][2] - a[0][0]*a[1][3]*a[2][2] - a[0][2]*a[1][0]*a[2][3] - a[0][3]*a[1][2]*a[2][0] ) / det;
  	M4D[2][0] = ( a[1][0]*a[2][1]*a[3][3] + a[1][1]*a[2][3]*a[3][0] + a[1][3]*a[2][0]*a[3][1] - a[1][0]*a[2][3]*a[3][1] - a[1][1]*a[2][0]*a[3][3] - a[1][3]*a[2][1]*a[3][0] ) / det;
  	M4D[2][1] = ( a[0][0]*a[2][3]*a[3][1] + a[0][1]*a[2][0]*a[3][3] + a[0][3]*a[2][1]*a[3][0] - a[0][0]*a[2][1]*a[3][3] - a[0][1]*a[2][3]*a[3][0] - a[0][3]*a[2][0]*a[3][1] ) / det;
  	M4D[2][2] = ( a[0][0]*a[1][1]*a[3][3] + a[0][1]*a[1][3]*a[3][0] + a[0][3]*a[1][0]*a[3][1] - a[0][0]*a[1][3]*a[3][1] - a[0][1]*a[1][0]*a[3][3] - a[0][3]*a[1][1]*a[3][0] ) / det;
  	M4D[2][3] = ( a[0][0]*a[1][3]*a[2][1] + a[0][1]*a[1][0]*a[2][3] + a[0][3]*a[1][1]*a[2][0] - a[0][0]*a[1][1]*a[2][3] - a[0][1]*a[1][3]*a[2][0] - a[0][3]*a[1][0]*a[2][1] ) / det;
  	M4D[3][0] = ( a[1][0]*a[2][2]*a[3][1] + a[1][1]*a[2][0]*a[3][2] + a[1][2]*a[2][1]*a[3][0] - a[1][0]*a[2][1]*a[3][2] - a[1][1]*a[2][2]*a[3][0] - a[1][2]*a[2][0]*a[3][1] ) / det;
  	M4D[3][1] = ( a[0][0]*a[2][1]*a[3][2] + a[0][1]*a[2][2]*a[3][0] + a[0][2]*a[2][0]*a[3][1] - a[0][0]*a[2][2]*a[3][1] - a[0][1]*a[2][0]*a[3][2] - a[0][2]*a[2][1]*a[3][0] ) / det;
  	M4D[3][2] = ( a[0][0]*a[1][2]*a[3][1] + a[0][1]*a[1][0]*a[3][2] + a[0][2]*a[1][1]*a[3][0] - a[0][0]*a[1][1]*a[3][2] - a[0][1]*a[1][2]*a[3][0] - a[0][2]*a[1][0]*a[3][1] ) / det;
  	M4D[3][3] = ( a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1] - a[0][0]*a[1][2]*a[2][1] - a[0][1]*a[1][0]*a[2][2] - a[0][2]*a[1][1]*a[2][0] ) / det;

  return det;
}

//逆行列を求める(掃出し法)整理しました。返り値はdetではない。返り値1:逆行列が存在する。返り値0:逆行列が存在しない。
int Inverse(double M[N_P_2][N_P_2])
{
  	int i, j, k, det_check=1;
  	double a[N_P_2][N_P_2], Tri[N_P_2][N_P_2];
  	double b[N_P_2][N_P_2]={0};
  	double check[N_P_2][N_P_2]={0};
  	double AKK, AIK;
  	double det = 1.0, buf;
  	int ErrorChecker=0;

  	/************行列式の計算。detとって逆行列の有無をチェック*************/
  	for(i=0; i<N_P_2; i++){
		for(j=0; j<N_P_2; j++){
	  		Tri[i][j] = M[i][j];
		}
  	}
  	//三角行列を作成
  	for(i=0;i<N_P_2;i++){
		for(j=0;j<N_P_2;j++){
	  		if(i<j){
				buf=Tri[j][i]/Tri[i][i];
				for(k=0;k<N_P_2;k++){
		  			Tri[j][k]-=Tri[i][k]*buf;
				}
	  		}
		}
  	}
  	//対角部分の積
  	for(i=0;i<N_P_2;i++){
		det*=Tri[i][i];
  	}
	if(fabs(det) < pow(10,-9)) det_check=0;


  	/*********************:本チャンの逆行列計算********************/
  	for(i=0; i<N_P_2; i++){
		for(j=0; j<N_P_2; j++){
	  		a[i][j] = M[i][j];;
	  		if(i==j) b[i][j] = 1.0;
		}
  	}

 	for(k=0; k<N_P_2; k++){
		AKK = a[k][k];
		for(j=0; j<N_P_2; j++){
	  		a[k][j] = a[k][j] / AKK;
	  		b[k][j] = b[k][j] / AKK;
		}
		for(i=0; i<N_P_2; i++){
	  		if (i != k){
	  			AIK = a[i][k];
	  			for(j=0; j<N_P_2; j++){
					a[i][j] = a[i][j] - (a[k][j] * AIK);
					b[i][j] = b[i][j] - (b[k][j] * AIK);
	  			}
			}
		}
 	}
  	//入力した行列と逆行列との積が単位行列になるかチェック
  	for(i=0; i<N_P_2; i++){
		for(j=0; j<N_P_2; j++){
	  		for(k=0; k<N_P_2; k++){
				check[i][j] += M[i][k] * b[k][j];
	  		}
		}
  	}

  	for(i=0;i<N_P_2;i++){
		for(j=0;j<N_P_2;j++){
			if(i==j)	{
				if(fabs(check[i][j]-1)>0.0000000001) printf("ERROR Failed to create InverseMatrix [%d][%d]!=1  (%f) \n",i,j,check[i][j] );
			}
			else{
				if(fabs(check[i][j])>0.0000000001) printf("ERROR Failed to create InverseMatrix [%d][%d]!=0  (%f)\n",i,j,check[i][j] );
			}
		}	
	}

  	for(i=0; i<N_P_2; i++){
		for(j=0; j<N_P_2; j++){
	  		M[i][j] = b[i][j];
		}
  	}
  	return det_check;

}

////////////////////////////////////////////////////////////////
//////////////////J_Integral////////////////////////////////////
////////////////////////////////////////////////////////////////

//1stP-Kを求める関数
void Make_1st_PK_stress(){

  	int i, j, k;
  	double K_D1[DIMENSION] = {1.0, 0.0, 0.0};
  	double K_D2[DIMENSION] = {0.0, 1.0, 0.0};
  	double K_D3[DIMENSION] = {0.0, 0.0, 1.0};
  	int e, N;

  	//変形勾配初期化
  	for( e = 0; e < NElements; e++ ){
		for( N = 0; N < N_IntegrationPoint; N++ ){
	  		for( i = 0; i < DIMENSION; i++ ){
				for( j = 0; j < DIMENSION; j++ ){
	  				Def_grad[e][N][i][j]=0.0;
				}
	  		}
		}
  	}

  	//変形勾配    I+∇u
 	for( e = 0; e < NElements; e++ ){
		for( N = 0; N < N_IntegrationPoint; N++ ){
	  		for( i = 0; i < DIMENSION; i++ ){
	  			Def_grad[e][N][i][0]=K_D1[i]+Disp_grad_1[e][N][i];
	  			Def_grad[e][N][i][1]=K_D2[i]+Disp_grad_2[e][N][i];
	  			Def_grad[e][N][i][2]=K_D3[i]+Disp_grad_3[e][N][i];
	  		}
		}
  	}
  	//ヤコビアンの計算
 	for( e = 0; e < NElements; e++ ){
   		for( N = 0; N < N_IntegrationPoint; N++ ){
	 		JJ[e][N]=Def_grad[e][N][0][0]*Def_grad[e][N][1][1]*Def_grad[e][N][2][2]
			 		+Def_grad[e][N][0][1]*Def_grad[e][N][1][2]*Def_grad[e][N][2][0]
			 		+Def_grad[e][N][0][2]*Def_grad[e][N][1][0]*Def_grad[e][N][2][1]
					-Def_grad[e][N][0][0]*Def_grad[e][N][1][2]*Def_grad[e][N][2][1]
					-Def_grad[e][N][0][1]*Def_grad[e][N][1][0]*Def_grad[e][N][2][2]
					-Def_grad[e][N][0][2]*Def_grad[e][N][1][1]*Def_grad[e][N][2][0];
  	 	}
 	}
  	//変形勾配の逆行列
 	Make_Def_grad_In();
 	//Output_Def_grad_Data();
  	//1stP-K
 	for( e = 0; e < NElements; e++ ){
		for( N = 0; N < N_IntegrationPoint; N++ ){
	  		for( i = 0; i < DIMENSION; i++ ){
				for( j = 0; j < DIMENSION; j++ ){
	  				if(deform==1){
						Pai[e][N][i][j]=Stress[e][N][i][j];
						//printf("Pai=%e\n",Pai[e][N][i][j]);
						//printf("Stress=%e\n",Stress[e][N][i][j]);
	  				}
	  				if(deform==2){  //1stPK=Jacobian * Stress * F^-1
						Pai[e][N][i][j]=JJ[e][N]*(Def_grad_In[e][N][i][0]*Stress[e][N][0][j]
					 					+Def_grad_In[e][N][i][1]*Stress[e][N][1][j]
					 					+Def_grad_In[e][N][i][2]*Stress[e][N][2][j]);
	  				}
				}
	  		}
		}
 	}
 //Output_Pai_Data();
}
//変形勾配の逆行列を求める
int Make_Def_grad_In()
{
  	int i, j, e, N;
  	double det;

  	for( e = 0; e < NElements; e++ ){
		for( N = 0; N < N_IntegrationPoint; N++ ){
	  		det=Def_grad[e][N][0][0]*Def_grad[e][N][1][1]*Def_grad[e][N][2][2]
				+Def_grad[e][N][0][1]*Def_grad[e][N][1][2]*Def_grad[e][N][2][0]
				+Def_grad[e][N][0][2]*Def_grad[e][N][1][0]*Def_grad[e][N][2][1]
				-Def_grad[e][N][0][0]*Def_grad[e][N][1][2]*Def_grad[e][N][2][1]
				-Def_grad[e][N][0][1]*Def_grad[e][N][1][0]*Def_grad[e][N][2][2]
				-Def_grad[e][N][0][2]*Def_grad[e][N][1][1]*Def_grad[e][N][2][0];

	 		Def_grad_In[e][N][0][0]=(Def_grad[e][N][1][1]*Def_grad[e][N][2][2]-Def_grad[e][N][1][2]*Def_grad[e][N][2][1])/det;
	 		Def_grad_In[e][N][0][1]=(Def_grad[e][N][2][1]*Def_grad[e][N][0][2]-Def_grad[e][N][0][1]*Def_grad[e][N][2][2])/det;
	  		Def_grad_In[e][N][0][2]=(Def_grad[e][N][0][1]*Def_grad[e][N][1][2]-Def_grad[e][N][1][1]*Def_grad[e][N][0][2])/det;
	  		Def_grad_In[e][N][1][0]=(Def_grad[e][N][2][0]*Def_grad[e][N][1][2]-Def_grad[e][N][1][0]*Def_grad[e][N][2][2])/det;
	  		Def_grad_In[e][N][1][1]=(Def_grad[e][N][0][0]*Def_grad[e][N][2][2]-Def_grad[e][N][2][0]*Def_grad[e][N][0][2])/det;
	  		Def_grad_In[e][N][1][2]=(Def_grad[e][N][1][0]*Def_grad[e][N][0][2]-Def_grad[e][N][0][0]*Def_grad[e][N][1][2])/det;
	  		Def_grad_In[e][N][2][0]=(Def_grad[e][N][1][0]*Def_grad[e][N][2][1]-Def_grad[e][N][2][0]*Def_grad[e][N][1][1])/det;
	  		Def_grad_In[e][N][2][1]=(Def_grad[e][N][2][0]*Def_grad[e][N][0][1]-Def_grad[e][N][0][0]*Def_grad[e][N][2][1])/det;
	  		Def_grad_In[e][N][2][2]=(Def_grad[e][N][0][0]*Def_grad[e][N][1][1]-Def_grad[e][N][1][0]*Def_grad[e][N][0][1])/det;
		}
  	}

  //for( e = 0; e < NElements; e++ ){
  //for( N = 0; N < N_IntegrationPoint; N++ ){
  //  for( i= 0; i< 3; i++ ){
  //    for( j= 0; j< 3; j++ ){
  //      I[e][N][i][j]=Def_grad_In[e][N][i][0]*Def_grad[e][N][0][j]
  //        +Def_grad_In[e][N][i][1]*Def_grad[e][N][1][j]
  //        +Def_grad_In[e][N][i][2]*Def_grad[e][N][2][j];
  //    }
  //  }
  //}
  //}
  Output_I();

  Output_Def_grad_In_Data();
}
/***********************************************************************
				要素のガウス点のglobal座標を求める。
***********************************************************************/
void make_Gxi_coord_Tetra_10()
{
	int e,N,i,j,k,l;
	double X[NO_NODES_ON_ELEMENT_T10][DIMENSION];
	double Gxi[4][3]={{0.13819660112501051518,0.58541019662496845446,0.13819660112501051518},
					 {0.13819660112501051518,0.13819660112501051518,0.13819660112501051518},
					 {0.13819660112501051518,0.13819660112501051518,0.58541019662496845446},
					 {0.58541019662496845446,0.13819660112501051518,0.13819660112501051518}};

	for( e = 0; e < NElements; e++ ){
		for( N = 0; N < N_IntegrationPoint; N++)
		{
			for( i = 0; i < DIMENSION; i++)
			{
				Gxi_coord[e][N][i] = 0.0;                    //初期化
			}
			for( i = 0; i < Number_of_Nodes_on_Elements; i++ )
			{
				for( j = 0; j < DIMENSION; j++ )
				{
					X[i][j] = Node_Coordinate[ ElementNodeId_s[e][i] ][j];  //要素ごとの節点IDの座標を格納
				}
			}
				for( i = 0; i < DIMENSION; i++ )
				{
					for(k = 0; k<  Number_of_Nodes_on_Elements; k++)
					{
						Gxi_coord[e][N][i] += N_Tetra_10( k, Gxi[N] ) * X[k][i];  //それぞれの要素のガウス点の座標をglobal座標系で表す
					}
				}
		}
	}
}
void make_Gxi_coord_Tetra_4()
{
	int e,N,i,j,k,l;
	double X[NO_NODES_ON_ELEMENT_T4][DIMENSION];
	double Gxi[1][3]={0.25,0.25,0.25};

	for( e = 0; e < NElements; e++ ){
		for( N = 0; N < N_IntegrationPoint; N++)
		{
			for( i = 0; i < DIMENSION; i++)
			{
				Gxi_coord[e][N][i] = 0.0;                    //初期化
			}
			for( i = 0; i < Number_of_Nodes_on_Elements; i++ )
			{
				for( j = 0; j < DIMENSION; j++ )
				{
					X[i][j] = Node_Coordinate[ ElementNodeId_s[e][i] ][j];  //要素ごとの節点IDの座標を格納
				}
			}
				for( i = 0; i < DIMENSION; i++ )
				{
					for(k = 0; k<  Number_of_Nodes_on_Elements; k++)
					{
						Gxi_coord[e][N][i] += N_Tetra_4( k, Gxi[N] ) * X[k][i];  //それぞれの要素のガウス点の座標をglobal座標系で表す
					}
				}
		}
	}
}
void make_Gxi_coord_Hexa_8()
{
	int e,N,i,j,k,l;
	double X[NO_NODES_ON_ELEMENT_H8][DIMENSION];
	double Gxi[8][3]=	{{-0.57735026918962584208,-0.57735026918962584208,-0.57735026918962584208},
					 { 0.57735026918962584208,-0.57735026918962584208,-0.57735026918962584208},
					 {-0.57735026918962584208, 0.57735026918962584208,-0.57735026918962584208},
					 { 0.57735026918962584208, 0.57735026918962584208,-0.57735026918962584208},
					 {-0.57735026918962584208,-0.57735026918962584208, 0.57735026918962584208},
					 { 0.57735026918962584208,-0.57735026918962584208, 0.57735026918962584208},
					 {-0.57735026918962584208, 0.57735026918962584208, 0.57735026918962584208},
					 { 0.57735026918962584208, 0.57735026918962584208, 0.57735026918962584208}};

	for( e = 0; e < NElements; e++ ){
		for( N = 0; N < N_IntegrationPoint; N++)
		{
			for( i = 0; i < DIMENSION; i++)
			{
				Gxi_coord[e][N][i] = 0.0;                    //初期化
			}
			for( i = 0; i < Number_of_Nodes_on_Elements; i++ )
			{
				for( j = 0; j < DIMENSION; j++ )
				{
					X[i][j] = Node_Coordinate[ ElementNodeId_s[e][i] ][j];  //要素ごとの節点IDの座標を格納
				}
			}		
			for( i = 0; i < DIMENSION; i++ )
			{
				for(k = 0; k<  Number_of_Nodes_on_Elements; k++)
				{
					Gxi_coord[e][N][i] += N_Hexa_8( k, Gxi[N] ) * X[k][i];  //それぞれの要素のガウス点の座標をglobal座標系で表す
				}
			}
		}
	}
}
void make_Gxi_coord_Hexa_20()
{
	int e,N,i,j,k,l;
	double X[NO_NODES_ON_ELEMENT_H20][DIMENSION];
	double Gxi[27][3]={{-0.77459666924148340428 	,-0.77459666924148340428 	,-0.77459666924148340428},
					  {  		0 	 		,-0.77459666924148340428 	,-0.77459666924148340428},
					  { 0.77459666924148340428 	,-0.77459666924148340428   	,-0.77459666924148340428},
					  {-0.77459666924148340428 	, 		 	0			,-0.77459666924148340428},
					  {  		0			,			0			,-0.77459666924148340428},
					  { 0.77459666924148340428 	, 			0			,-0.77459666924148340428},
					  {-0.77459666924148340428 	, 0.77459666924148340428 	,-0.77459666924148340428},
					  { 		 	0 			, 0.77459666924148340428 	,-0.77459666924148340428},
					  { 0.77459666924148340428 	, 0.77459666924148340428 	,-0.77459666924148340428},
					  {-0.77459666924148340428 	,-0.77459666924148340428 	,		 0			},
					  {  		0 			,-0.77459666924148340428 	, 		 0			},
					  { 0.77459666924148340428 	,-0.77459666924148340428   	, 		 0			},
					  {-0.77459666924148340428 	,  			0			,  		 0			},
					  {  		0			, 			0			,  		 0   		},
					  { 0.77459666924148340428 	,  			0			,  		 0			},
					  {-0.77459666924148340428 	, 0.77459666924148340428 	,  		 0			},
					  {  		0		 	, 0.77459666924148340428 	,  		 0			},
					  { 0.77459666924148340428 	, 0.77459666924148340428 	,  		 0			},
					  {-0.77459666924148340428 	,-0.77459666924148340428 	, 0.77459666924148340428},
					  {  		0 			,-0.77459666924148340428 	, 0.77459666924148340428},
					  { 0.77459666924148340428 	,-0.77459666924148340428   	, 0.77459666924148340428},
					  {-0.77459666924148340428 	,  			0			, 0.77459666924148340428},
					  {  		0			, 			0			, 0.77459666924148340428},
					  { 0.77459666924148340428 	, 			0			, 0.77459666924148340428},
					  {-0.77459666924148340428 	, 0.77459666924148340428 	, 0.77459666924148340428},
					  { 		 	0 			, 0.77459666924148340428 	, 0.77459666924148340428},
					  { 0.77459666924148340428 	, 0.77459666924148340428 	, 0.77459666924148340428}};

	for( e = 0; e < NElements; e++ ){
		for( N = 0; N < N_IntegrationPoint; N++)
		{
			for( i = 0; i < DIMENSION; i++)
			{
				Gxi_coord[e][N][i] = 0.0;                    //初期化
			}
			for( i = 0; i < Number_of_Nodes_on_Elements; i++ )
			{
				for( j = 0; j < DIMENSION; j++ )
				{
					X[i][j] = Node_Coordinate[ ElementNodeId_s[e][i] ][j];  //要素ごとの節点IDの座標を格納
				}
			}
				for( i = 0; i < DIMENSION; i++ )
				{
					for(k = 0; k<  Number_of_Nodes_on_Elements; k++)
					{
						Gxi_coord[e][N][i] += N_Hexa_20( k, Gxi[N] ) * X[k][i];  //それぞれの要素のガウス点の座標をglobal座標系で表す
					}
				}
		}
	}
}
//ひずみエネルギー密度直接or体積変化を考慮したものをJWに代入する。
int Make_JW()
{
  	int i,k;
  	//ひずみエネルギー密度に体積変化を考慮した。
  	if(deform==1){
  		for( i = 0; i < NElements; i++ ){
  			for( k = 0; k < N_IntegrationPoint; k++ ){
	  			JW[i][k]= W[i][k];
			}
  		}
  	}
  	else if(deform==2){
  		for( i = 0; i < NElements; i++ ){
  			for( k = 0; k < N_IntegrationPoint; k++ ){
	  			JW[i][k]= JJ[i][k] * W[i][k];
			}
  		}
  	}
}
//エネルギーモーメンタムテンソルを求める関数
void Make_EMT(){

  	int i, j, k;
  	double K_D1[DIMENSION] = {1.0, 0.0, 0.0};
  	double K_D2[DIMENSION] = {0.0, 1.0, 0.0};
  	double K_D3[DIMENSION] = {0.0, 0.0, 1.0};

  	//謎の初期化
  	for( k = 0; k < N_IntegrationPoint; k++ ){
		for( i = 0; i < NElements; i++ ){
	  		for( j = 0; j < DIMENSION; j++ ){
				W_1[i][k][j] = 0.0;
				W_2[i][k][j] = 0.0;
				W_3[i][k][j] = 0.0;
	  		}
		}
  	}

  //体積変化を考慮したひずみエネルギー密度とクロネッカーデルタの積をとる。
  	for( k = 0; k < N_IntegrationPoint; k++ ){
		for( i = 0; i < NElements; i++ ){
	  		for( j = 0; j < DIMENSION; j++ ){
				W_K_D1[i][k][j] = JW[i][k] * K_D1[j];
				W_K_D2[i][k][j] = JW[i][k] * K_D2[j];
				W_K_D3[i][k][j] = JW[i][k] * K_D3[j];
	  		}
		}
  	}

  	for( k = 0; k < N_IntegrationPoint; k++ ){
		for( i = 0; i < NElements; i++ ){
	
			W_1[i][k][0] += Pai[i][k][0][0] * Disp_grad_1[i][k][0]
	  						+ Pai[i][k][0][1] * Disp_grad_1[i][k][1] + Pai[i][k][0][2] * Disp_grad_1[i][k][2];
	
			W_1[i][k][1] += Pai[i][k][1][0] * Disp_grad_1[i][k][0]
	  						+ Pai[i][k][1][1] * Disp_grad_1[i][k][1] + Pai[i][k][1][2] * Disp_grad_1[i][k][2];
	
			W_1[i][k][2] += Pai[i][k][2][0] * Disp_grad_1[i][k][0]
	  						+ Pai[i][k][2][1] * Disp_grad_1[i][k][1] + Pai[i][k][2][2] * Disp_grad_1[i][k][2];
	
			W_2[i][k][0] += Pai[i][k][0][0] * Disp_grad_2[i][k][0]
	  						+ Pai[i][k][0][1] * Disp_grad_2[i][k][1] + Pai[i][k][0][2] * Disp_grad_2[i][k][2];
	
			W_2[i][k][1] += Pai[i][k][1][0] * Disp_grad_2[i][k][0]
	  						+ Pai[i][k][1][1] * Disp_grad_2[i][k][1] + Pai[i][k][1][2] * Disp_grad_2[i][k][2];
	
			W_2[i][k][2] += Pai[i][k][2][0] * Disp_grad_2[i][k][0]
	  						+ Pai[i][k][2][1] * Disp_grad_2[i][k][1] + Pai[i][k][2][2] * Disp_grad_2[i][k][2];
		
			W_3[i][k][0] += Pai[i][k][0][0] * Disp_grad_3[i][k][0]
	  						+ Pai[i][k][0][1] * Disp_grad_3[i][k][1] + Pai[i][k][0][2] * Disp_grad_3[i][k][2];
	
			W_3[i][k][1] += Pai[i][k][1][0] * Disp_grad_3[i][k][0]
	  						+ Pai[i][k][1][1] * Disp_grad_3[i][k][1] + Pai[i][k][1][2] * Disp_grad_3[i][k][2];
	
			W_3[i][k][2] += Pai[i][k][2][0] * Disp_grad_3[i][k][0]
	  						+ Pai[i][k][2][1] * Disp_grad_3[i][k][1] + Pai[i][k][2][2] * Disp_grad_3[i][k][2];
		}
  	}
}

/***********************************************************************
全要素での変位勾配作成
***********************************************************************/
void Make_Disp_grad_2_Tetra_10(){
	static double U[KIEL_SIZE_T10];
	static double B_1[DIMENSION][KIEL_SIZE_T10], B_2[DIMENSION][KIEL_SIZE_T10], B_3[DIMENSION][KIEL_SIZE_T10];
	static double X[NO_NODES_ON_ELEMENT_T10][DIMENSION],J;
	int N,e,i,j;
	double Gxi[4][3]={{0.13819660112501051518,0.58541019662496845446,0.13819660112501051518},
					 {0.13819660112501051518,0.13819660112501051518,0.13819660112501051518},
					 {0.13819660112501051518,0.13819660112501051518,0.58541019662496845446},
					 {0.58541019662496845446,0.13819660112501051518,0.13819660112501051518}};
	int k,l;

	if(DispMode==0){
		for( e = 0; e < NElements; e++ ){
			for( i = 0; i < Number_of_Nodes_on_Elements; i++ ){
				for( j = 0; j < DIMENSION; j++ ){
					U[ i*DIMENSION +j ] = Displacement[ElementNodeId_s[e][i]][j];
					X[i][j] = Node_Coordinate[ ElementNodeId_s[e][i] ][j];      //初期配置
					//printf("%lf ", X[i][j] );
				}
			}

			for( N = 0; N < N_IntegrationPoint; N++ ){
				//printf("e=%d\n",e );
		  		Make_Gradient_Matrix_Tetra_10(B_1, B_2, B_3, Gxi[N], X ,&J);
		  		//printf("e=%d\n",e );
				for( i = 0; i < DIMENSION; i++ ){
					for( j = 0; j < KIEL_SIZE_T10; j++ ){
				  		Disp_grad_1[e][N][i] += B_1[i][j] * U[j]; //Displacementの値の∂N/∂x成分を計算
				  		Disp_grad_2[e][N][i] += B_2[i][j] * U[j]; //Displacementの値の∂N/∂y成分を計算
				  		Disp_grad_3[e][N][i] += B_3[i][j] * U[j]; //Displacementの値の∂N/∂z成分を計算
				  		//printf("%lf %lf %lf |%lf\n",B_1[i][j],B_2[i][j],B_3[i][j],U[j] );
					}
					//printf("DispGlad ele %d _ %d _%d |%lf %lf %lf\n",e,N,i,Disp_grad_1[e][N][i],Disp_grad_2[e][N][i],Disp_grad_3[e][N][i] );
				}
			}
		}
	}
	else if(DispMode==1){
		for( e = 0; e < NElements; e++ ){
			for( i = 0; i < Number_of_Nodes_on_Elements; i++ ){
				for( j = 0; j < DIMENSION; j++ ){
					U[ i*DIMENSION +j ] = Displacement[ ElementNodeId_s[e][i]][j];
					X[i][j] = Node_Coordinate[ ElementNodeId_s[e][i] ][j]
									+Displacement[ ElementNodeId_s[e][i]][j];  //現在配置
				}
			}
			for( N = 0; N < N_IntegrationPoint; N++ ){
		  		Make_Gradient_Matrix_Tetra_10(B_1, B_2, B_3, Gxi[N], X ,&J);
				for( i = 0; i < DIMENSION; i++ ){
					for( j = 0; j < KIEL_SIZE_T10; j++ ){
				  		Disp_grad_1[e][N][i] += B_1[i][j] * U[j]; //Displacementの値の∂N/∂x成分を計算
				  		Disp_grad_2[e][N][i] += B_2[i][j] * U[j]; //Displacementの値の∂N/∂y成分を計算
				  		Disp_grad_3[e][N][i] += B_3[i][j] * U[j]; //Displacementの値の∂N/∂z成分を計算
					}
				}
			}
		}
	}

	Output_Disp_grad_Data();
}
void Make_Disp_grad_2_Tetra_4(){
	static double U[KIEL_SIZE_T4];
	static double B_1[DIMENSION][KIEL_SIZE_T4], B_2[DIMENSION][KIEL_SIZE_T4], B_3[DIMENSION][KIEL_SIZE_T4];
	static double X[NO_NODES_ON_ELEMENT_T4][DIMENSION],J;
	int N,e,i,j;
	double Gxi[1][3]={0.25,0.25,0.25};

	if(DispMode==0){
		for( e = 0; e < NElements; e++ ){
			for( i = 0; i < Number_of_Nodes_on_Elements; i++ ){
				for( j = 0; j < DIMENSION; j++ ){
					U[ i*DIMENSION +j ] = Displacement[ElementNodeId_s[e][i]][j];
					X[i][j] = Node_Coordinate[ ElementNodeId_s[e][i] ][j];      //初期配置
				}
			}
			for( N = 0; N < N_IntegrationPoint; N++ ){
		  		Make_Gradient_Matrix_Tetra_4(B_1, B_2, B_3, Gxi[N], X ,&J );
				for( i = 0; i < DIMENSION; i++ ){
					for( j = 0; j < KIEL_SIZE_T4; j++ ){
				  		Disp_grad_1[e][N][i] += B_1[i][j] * U[j]; //Displacementの値の∂N/∂x成分を計算
				  		Disp_grad_2[e][N][i] += B_2[i][j] * U[j]; //Displacementの値の∂N/∂y成分を計算
				  		Disp_grad_3[e][N][i] += B_3[i][j] * U[j]; //Displacementの値の∂N/∂z成分を計算
					}
				}
			}
		}
	}
	else if(DispMode==1){
		for( e = 0; e < NElements; e++ ){
			for( i = 0; i < Number_of_Nodes_on_Elements; i++ ){
				for( j = 0; j < DIMENSION; j++ ){
					U[ i*DIMENSION +j ] = Displacement[ ElementNodeId_s[e][i]][j];
					X[i][j] = Node_Coordinate[ ElementNodeId_s[e][i] ][j]
									+Displacement[ ElementNodeId_s[e][i]][j];  //現在配置
				}
			}
			for( N = 0; N < N_IntegrationPoint; N++ ){
		  		Make_Gradient_Matrix_Tetra_4(B_1, B_2, B_3, Gxi[N], X ,&J );
				for( i = 0; i < DIMENSION; i++ ){
					for( j = 0; j < KIEL_SIZE_T4; j++ ){
				  		Disp_grad_1[e][N][i] += B_1[i][j] * U[j]; //Displacementの値の∂N/∂x成分を計算
				  		Disp_grad_2[e][N][i] += B_2[i][j] * U[j]; //Displacementの値の∂N/∂y成分を計算
				  		Disp_grad_3[e][N][i] += B_3[i][j] * U[j]; //Displacementの値の∂N/∂z成分を計算
					}
				}
			}
		}
	}
	//Output_Disp_grad_Data();
}
void Make_Disp_grad_2_Hexa_8(){
	static double U[KIEL_SIZE_H8];
	static double B_1[DIMENSION][KIEL_SIZE_H8], B_2[DIMENSION][KIEL_SIZE_H8], B_3[DIMENSION][KIEL_SIZE_H8];
	static double X[NO_NODES_ON_ELEMENT_H8][DIMENSION],J;
	int N,e,i,j;
	double Gxi[8][3]=	{{-0.57735026918962584208,-0.57735026918962584208,-0.57735026918962584208},
					 { 0.57735026918962584208,-0.57735026918962584208,-0.57735026918962584208},
					 {-0.57735026918962584208, 0.57735026918962584208,-0.57735026918962584208},
					 { 0.57735026918962584208, 0.57735026918962584208,-0.57735026918962584208},
					 {-0.57735026918962584208,-0.57735026918962584208, 0.57735026918962584208},
					 { 0.57735026918962584208,-0.57735026918962584208, 0.57735026918962584208},
					 {-0.57735026918962584208, 0.57735026918962584208, 0.57735026918962584208},
					 { 0.57735026918962584208, 0.57735026918962584208, 0.57735026918962584208}};

	if(DispMode==0){
		for( e = 0; e < NElements; e++ ){
			for( i = 0; i < Number_of_Nodes_on_Elements; i++ ){
				for( j = 0; j < DIMENSION; j++ ){
					U[ i*DIMENSION +j ] = Displacement[ElementNodeId_s[e][i]][j];
					X[i][j] = Node_Coordinate[ ElementNodeId_s[e][i] ][j];      //初期配置
				}
			}
			for( N = 0; N < N_IntegrationPoint; N++ ){
		  		Make_Gradient_Matrix_Hexa_8(B_1, B_2, B_3, Gxi[N], X ,&J );
				for( i = 0; i < DIMENSION; i++ ){
					for( j = 0; j < KIEL_SIZE_H8; j++ ){
				  		Disp_grad_1[e][N][i] += B_1[i][j] * U[j]; //Displacementの値の∂N/∂x成分を計算
				  		Disp_grad_2[e][N][i] += B_2[i][j] * U[j]; //Displacementの値の∂N/∂y成分を計算
				  		Disp_grad_3[e][N][i] += B_3[i][j] * U[j]; //Displacementの値の∂N/∂z成分を計算
					}
				}
			}
		}
		/*for( e = 0; e < NElements; e++ ){//kesuzooooooooooooooooo
			for( N = 0; N < N_IntegrationPoint; N++ ){
					printf("Ele= %d iGP= %d | %lf %lf %lf | %lf %lf %lf | %lf %lf %lf |\n",e,N,Disp_grad_1[e][N][0],
					Disp_grad_1[e][N][1],Disp_grad_1[e][N][2],Disp_grad_2[e][N][0],Disp_grad_2[e][N][1],Disp_grad_2[e][N][2],
					Disp_grad_3[e][N][0],Disp_grad_3[e][N][1],Disp_grad_3[e][N][2]);
			}
		}*/


	}
	else if(DispMode==1){
		for( e = 0; e < NElements; e++ ){
			for( i = 0; i < Number_of_Nodes_on_Elements; i++ ){
				for( j = 0; j < DIMENSION; j++ ){
					U[ i*DIMENSION +j ] = Displacement[ ElementNodeId_s[e][i]][j];
					X[i][j] = Node_Coordinate[ ElementNodeId_s[e][i] ][j]
									+Displacement[ ElementNodeId_s[e][i]][j];  //現在配置
				}
			}
			for( N = 0; N < N_IntegrationPoint; N++ ){
		  		Make_Gradient_Matrix_Hexa_8(B_1, B_2, B_3, Gxi[N], X ,&J );
				for( i = 0; i < DIMENSION; i++ ){
					for( j = 0; j < KIEL_SIZE_H8; j++ ){
				  		Disp_grad_1[e][N][i] += B_1[i][j] * U[j]; //Displacementの値の∂N/∂x成分を計算
				  		Disp_grad_2[e][N][i] += B_2[i][j] * U[j]; //Displacementの値の∂N/∂y成分を計算
				  		Disp_grad_3[e][N][i] += B_3[i][j] * U[j]; //Displacementの値の∂N/∂z成分を計算
					}
				}
			}
		}
	}
	Output_Disp_grad_Data();
}
void Make_Disp_grad_2_Hexa_20(){
	static double U[KIEL_SIZE_H20];
	static double B_1[DIMENSION][KIEL_SIZE_H20], B_2[DIMENSION][KIEL_SIZE_H20], B_3[DIMENSION][KIEL_SIZE_H20];
	static double X[NO_NODES_ON_ELEMENT_H20][DIMENSION],J;
	int N,e,i,j;
	double Gxi[27][3]={{-0.77459666924148340428 	,-0.77459666924148340428 	,-0.77459666924148340428},
					  {  		0 	 		,-0.77459666924148340428 	,-0.77459666924148340428},
					  { 0.77459666924148340428 	,-0.77459666924148340428   	,-0.77459666924148340428},
					  {-0.77459666924148340428 	, 		 	0			,-0.77459666924148340428},
					  {  		0			,			0			,-0.77459666924148340428},
					  { 0.77459666924148340428 	, 			0			,-0.77459666924148340428},
					  {-0.77459666924148340428 	, 0.77459666924148340428 	,-0.77459666924148340428},
					  { 		 	0 			, 0.77459666924148340428 	,-0.77459666924148340428},
					  { 0.77459666924148340428 	, 0.77459666924148340428 	,-0.77459666924148340428},
					  {-0.77459666924148340428 	,-0.77459666924148340428 	,		 0			},
					  {  		0 			,-0.77459666924148340428 	, 		 0			},
					  { 0.77459666924148340428 	,-0.77459666924148340428   	, 		 0			},
					  {-0.77459666924148340428 	,  			0			,  		 0			},
					  {  		0			, 			0			,  		 0   		},
					  { 0.77459666924148340428 	,  			0			,  		 0			},
					  {-0.77459666924148340428 	, 0.77459666924148340428 	,  		 0			},
					  {  		0		 	, 0.77459666924148340428 	,  		 0			},
					  { 0.77459666924148340428 	, 0.77459666924148340428 	,  		 0			},
					  {-0.77459666924148340428 	,-0.77459666924148340428 	, 0.77459666924148340428},
					  {  		0 			,-0.77459666924148340428 	, 0.77459666924148340428},
					  { 0.77459666924148340428 	,-0.77459666924148340428   	, 0.77459666924148340428},
					  {-0.77459666924148340428 	,  			0			, 0.77459666924148340428},
					  {  		0			, 			0			, 0.77459666924148340428},
					  { 0.77459666924148340428 	, 			0			, 0.77459666924148340428},
					  {-0.77459666924148340428 	, 0.77459666924148340428 	, 0.77459666924148340428},
					  { 		 	0 			, 0.77459666924148340428 	, 0.77459666924148340428},
					  { 0.77459666924148340428 	, 0.77459666924148340428 	, 0.77459666924148340428}};

	if(DispMode==0){
		for( e = 0; e < NElements; e++ ){
			for( i = 0; i < Number_of_Nodes_on_Elements; i++ ){
				for( j = 0; j < DIMENSION; j++ ){
					U[ i*DIMENSION +j ] = Displacement[ElementNodeId_s[e][i]][j];
					X[i][j] = Node_Coordinate[ ElementNodeId_s[e][i] ][j];      //初期配置
				}
			}
			for( N = 0; N < N_IntegrationPoint; N++ ){
		  		Make_Gradient_Matrix_Hexa_20(B_1, B_2, B_3, Gxi[N], X ,&J );
				for( i = 0; i < DIMENSION; i++ ){
					for( j = 0; j < KIEL_SIZE_H20; j++ ){
				  		Disp_grad_1[e][N][i] += B_1[i][j] * U[j]; //Displacementの値の∂N/∂x成分を計算
				  		Disp_grad_2[e][N][i] += B_2[i][j] * U[j]; //Displacementの値の∂N/∂y成分を計算
				  		Disp_grad_3[e][N][i] += B_3[i][j] * U[j]; //Displacementの値の∂N/∂z成分を計算
					}
				}
			}
		}
	}
	else if(DispMode==1){
		for( e = 0; e < NElements; e++ ){
			for( i = 0; i < Number_of_Nodes_on_Elements; i++ ){
				for( j = 0; j < DIMENSION; j++ ){
					U[ i*DIMENSION +j ] = Displacement[ ElementNodeId_s[e][i]][j];
					X[i][j] = Node_Coordinate[ ElementNodeId_s[e][i] ][j]
									+Displacement[ ElementNodeId_s[e][i]][j];  //現在配置
				}
			}
			for( N = 0; N < N_IntegrationPoint; N++ ){
		  		Make_Gradient_Matrix_Hexa_20(B_1, B_2, B_3, Gxi[N], X ,&J );
				for( i = 0; i < DIMENSION; i++ ){
					for( j = 0; j < KIEL_SIZE_H20; j++ ){
				  		Disp_grad_1[e][N][i] += B_1[i][j] * U[j]; //Displacementの値の∂N/∂x成分を計算
				  		Disp_grad_2[e][N][i] += B_2[i][j] * U[j]; //Displacementの値の∂N/∂y成分を計算
				  		Disp_grad_3[e][N][i] += B_3[i][j] * U[j]; //Displacementの値の∂N/∂z成分を計算
					}
				}
			}
		}
	}
	//Output_Disp_grad_Data();
}
/***********************************************************************
勾配を求めるための行列作成。※もともとBマトリクス作成ツールを改造したものらしい。
***********************************************************************/

int Make_Gradient_Matrix_Tetra_10(double B_1[DIMENSION][KIEL_SIZE_T10], double B_2[DIMENSION][KIEL_SIZE_T10], double B_3[DIMENSION][KIEL_SIZE_T10],
				 double Local_coord[DIMENSION], double X[NO_NODES_ON_ELEMENT_T10][DIMENSION], double *J )
{

  	static double a[DIMENSION][DIMENSION], b[DIMENSION][NO_NODES_ON_ELEMENT_T10];
	int i,j,k;

	Jacobian_Tetra_10( a, Local_coord, X );

	*J = InverseMatrix_3D( a )/6.0;
		//printf("%lf\n",*J);

	if( *J <= 0 ){
		return -999;
	}


	/*ヤコビアンの逆行列[3×3]×形状関数の微分[3×10(要素にある節点数)]*/
	for( i = 0; i < DIMENSION; i++ ){
		for( j = 0; j < Number_of_Nodes_on_Elements; j++ ){
			b[i][j] = 0.0;
			for( k = 0; k < DIMENSION; k++ ){
				b[i][j] += a[k][i] * dN_Tetra_10( j, Local_coord, k);
			}
			//printf("b[%d][%d]=%lf ",i,j,b[i][j] );
		}
		//printf("\n");
	}
	//B1=∂N/∂xの対角行列がN0~N10まで並んでる。B2は∂N/∂y...
	for( i = 0; i < Number_of_Nodes_on_Elements; i++ ){
	  	B_1[0][3*i] = b[0][i];    B_1[0][3*i+1] = 0.0;        B_1[0][3*i+2] = 0.0;
	  	B_1[1][3*i] = 0.0;        B_1[1][3*i+1] = b[0][i];    B_1[1][3*i+2] = 0.0;
	  	B_1[2][3*i] = 0.0;        B_1[2][3*i+1] = 0.0;        B_1[2][3*i+2] = b[0][i];

		B_2[0][3*i] = b[1][i];    B_2[0][3*i+1] = 0.0;        B_2[0][3*i+2] = 0.0;
	  	B_2[1][3*i] = 0.0;        B_2[1][3*i+1] = b[1][i];    B_2[1][3*i+2] = 0.0;
	  	B_2[2][3*i] = 0.0;        B_2[2][3*i+1] = 0.0;        B_2[2][3*i+2] = b[1][i];
	
	  	B_3[0][3*i] = b[2][i];    B_3[0][3*i+1] = 0.0;        B_3[0][3*i+2] = 0.0;
	  	B_3[1][3*i] = 0.0;        B_3[1][3*i+1] = b[2][i];    B_3[1][3*i+2] = 0.0;
	  	B_3[2][3*i] = 0.0;        B_3[2][3*i+1] = 0.0;        B_3[2][3*i+2] = b[2][i];		
	}

	return 0;
}

int Make_Gradient_Matrix_Tetra_4(double B_1[DIMENSION][KIEL_SIZE_T4], double B_2[DIMENSION][KIEL_SIZE_T4], double B_3[DIMENSION][KIEL_SIZE_T4],
				 double Local_coord[DIMENSION], double X[NO_NODES_ON_ELEMENT_T4][DIMENSION], double *J )
{
  	static double a[DIMENSION][DIMENSION], b[DIMENSION][NO_NODES_ON_ELEMENT_T4];
	int i,j,k;

	Jacobian_Tetra_4( a, Local_coord, X );
	*J = InverseMatrix_3D( a )/6.0;
	if( *J <= 0 )return -999;

	/*ヤコビアンの逆行列[3×3]×形状関数の微分[3×10(要素にある節点数)]*/
	for( i = 0; i < DIMENSION; i++ ){
		for( j = 0; j < Number_of_Nodes_on_Elements; j++ ){
			b[i][j] = 0.0;
			for( k = 0; k < DIMENSION; k++ ){
				b[i][j] += a[k][i] * dN_Tetra_4( j, Local_coord, k);
			}
		}
	}
	//B1=∂N/∂xの対角行列がN0~N10まで並んでる。B2は∂N/∂y...
	for( i = 0; i < Number_of_Nodes_on_Elements; i++ ){
	  	B_1[0][3*i] = b[0][i];    B_1[0][3*i+1] = 0.0;        B_1[0][3*i+2] = 0.0;
	  	B_1[1][3*i] = 0.0;        B_1[1][3*i+1] = b[0][i];    B_1[1][3*i+2] = 0.0;
	  	B_1[2][3*i] = 0.0;        B_1[2][3*i+1] = 0.0;        B_1[2][3*i+2] = b[0][i];

		B_2[0][3*i] = b[1][i];    B_2[0][3*i+1] = 0.0;        B_2[0][3*i+2] = 0.0;
	  	B_2[1][3*i] = 0.0;        B_2[1][3*i+1] = b[1][i];    B_2[1][3*i+2] = 0.0;
	  	B_2[2][3*i] = 0.0;        B_2[2][3*i+1] = 0.0;        B_2[2][3*i+2] = b[1][i];
	
	  	B_3[0][3*i] = b[2][i];    B_3[0][3*i+1] = 0.0;        B_3[0][3*i+2] = 0.0;
	  	B_3[1][3*i] = 0.0;        B_3[1][3*i+1] = b[2][i];    B_3[1][3*i+2] = 0.0;
	  	B_3[2][3*i] = 0.0;        B_3[2][3*i+1] = 0.0;        B_3[2][3*i+2] = b[2][i];	
	}

	return 0;
}

int Make_Gradient_Matrix_Hexa_8(double B_1[DIMENSION][KIEL_SIZE_H8], double B_2[DIMENSION][KIEL_SIZE_H8], double B_3[DIMENSION][KIEL_SIZE_H8],
				 double Local_coord[DIMENSION], double X[NO_NODES_ON_ELEMENT_H8][DIMENSION], double *J )
{

  	static double a[DIMENSION][DIMENSION], b[DIMENSION][NO_NODES_ON_ELEMENT_H8];
	int i,j,k;

	Jacobian_Hexa_8( a, Local_coord, X );
	*J = InverseMatrix_3D( a );
	if( *J <= 0 )return -999;

	/*ヤコビアンの逆行列[3×3]×形状関数の微分[3×10(要素にある節点数)]*/
	for( i = 0; i < DIMENSION; i++ ){
		for( j = 0; j < Number_of_Nodes_on_Elements; j++ ){
			b[i][j] = 0.0;
			for( k = 0; k < DIMENSION; k++ ){
				b[i][j] += a[k][i] * dN_Hexa_8( j, Local_coord, k);
			}
		}
	}
	//B1=∂N/∂xの対角行列がN0~N10まで並んでる。B2は∂N/∂y...
	for( i = 0; i < Number_of_Nodes_on_Elements; i++ ){
	  	B_1[0][3*i] = b[0][i];    B_1[0][3*i+1] = 0.0;        B_1[0][3*i+2] = 0.0;
	  	B_1[1][3*i] = 0.0;        B_1[1][3*i+1] = b[0][i];    B_1[1][3*i+2] = 0.0;
	  	B_1[2][3*i] = 0.0;        B_1[2][3*i+1] = 0.0;        B_1[2][3*i+2] = b[0][i];

		B_2[0][3*i] = b[1][i];    B_2[0][3*i+1] = 0.0;        B_2[0][3*i+2] = 0.0;
	  	B_2[1][3*i] = 0.0;        B_2[1][3*i+1] = b[1][i];    B_2[1][3*i+2] = 0.0;
	  	B_2[2][3*i] = 0.0;        B_2[2][3*i+1] = 0.0;        B_2[2][3*i+2] = b[1][i];
	
	  	B_3[0][3*i] = b[2][i];    B_3[0][3*i+1] = 0.0;        B_3[0][3*i+2] = 0.0;
	  	B_3[1][3*i] = 0.0;        B_3[1][3*i+1] = b[2][i];    B_3[1][3*i+2] = 0.0;
	  	B_3[2][3*i] = 0.0;        B_3[2][3*i+1] = 0.0;        B_3[2][3*i+2] = b[2][i];	
	}
	//printf("B_1: %lf, B_2: %lf, B_3: %lf\n", B_1[0][0], B_2[0][0], B_3[0][0]);
	return 0;
}
int Make_Gradient_Matrix_Hexa_20(double B_1[DIMENSION][KIEL_SIZE_H20], double B_2[DIMENSION][KIEL_SIZE_H20], double B_3[DIMENSION][KIEL_SIZE_H20],
				 double Local_coord[DIMENSION], double X[NO_NODES_ON_ELEMENT_H20][DIMENSION], double *J )
{
  	static double a[DIMENSION][DIMENSION], b[DIMENSION][NO_NODES_ON_ELEMENT_H20];
	int i,j,k;

	Jacobian_Hexa_20( a, Local_coord, X );
	*J = InverseMatrix_3D( a );
	if( *J <= 0 )return -999;

	/*ヤコビアンの逆行列[3×3]×形状関数の微分[3×10(要素にある節点数)]*/
	for( i = 0; i < DIMENSION; i++ ){
		for( j = 0; j < Number_of_Nodes_on_Elements; j++ ){
			b[i][j] = 0.0;
			for( k = 0; k < DIMENSION; k++ ){
				b[i][j] += a[k][i] * dN_Hexa_20( j, Local_coord, k);
			}
		}
	}
	//B1=∂N/∂xの対角行列がN0~N10まで並んでる。B2は∂N/∂y...
	for( i = 0; i < Number_of_Nodes_on_Elements; i++ ){
	  	B_1[0][3*i] = b[0][i];    B_1[0][3*i+1] = 0.0;        B_1[0][3*i+2] = 0.0;
	  	B_1[1][3*i] = 0.0;        B_1[1][3*i+1] = b[0][i];    B_1[1][3*i+2] = 0.0;
	  	B_1[2][3*i] = 0.0;        B_1[2][3*i+1] = 0.0;        B_1[2][3*i+2] = b[0][i];

		B_2[0][3*i] = b[1][i];    B_2[0][3*i+1] = 0.0;        B_2[0][3*i+2] = 0.0;
	  	B_2[1][3*i] = 0.0;        B_2[1][3*i+1] = b[1][i];    B_2[1][3*i+2] = 0.0;
	  	B_2[2][3*i] = 0.0;        B_2[2][3*i+1] = 0.0;        B_2[2][3*i+2] = b[1][i];
	
	  	B_3[0][3*i] = b[2][i];    B_3[0][3*i+1] = 0.0;        B_3[0][3*i+2] = 0.0;
	  	B_3[1][3*i] = 0.0;        B_3[1][3*i+1] = b[2][i];    B_3[1][3*i+2] = 0.0;
	  	B_3[2][3*i] = 0.0;        B_3[2][3*i+1] = 0.0;        B_3[2][3*i+2] = b[2][i];	
	}

	return 0;
}

/***********************************************************************
						最小二乗近似
***********************************************************************/
//最小二乗近似(線形多項式近似) (P(応力等の数値) = a + bx + cy + dz)で近似している。
int LSM(int ip, int N_SP, int data, int N_COMPONENT)
{
  	//変数の設定
  	double a[N_COMPONENT][N_P]; //近似多項式の係数(これを求める)
  	double P[N_SP][N_P]; //基底
  	double A[N_P][N_P]; //近似多項式をa[i]についての偏微分してできる連立一次方程式におけるaの係数行列
  	double A_Inv[N_P][N_P]; //行列Aの逆行列
  	double b[N_COMPONENT][N_P]; //近似多項式をa[i]についての偏微分してできる連立一次方程式における右辺の項
  	double check[N_P][N_P];
  	int i, j, k, l;
  	
  	//配列の初期化
  	for(j=0; j<N_SP; j++){
		for(k=0; k<N_P; k++){
	  		P[j][k] = 0.0;
		}
  	}
  	for(i=0; i<N_COMPONENT; i++){
		for(k=0; k<N_P; k++){
	  		b[i][k] = 0.0;
		}
  	}
  	//基底ベクトルP(線形近似_n次元){1,x,y,z}がサンプリング点数だけある形の行列。
  	for(j=0; j<N_SP; j++){
		P[j][0] = 1.0;
		for(k=1; k<N_P; k++){
	  		P[j][k] = Coord_of_SP[j][k-1]-Coord_of_PatchCenter[ip][k-1];
		}
  	}

  	//[Pの行列][ガウス点の持つ値の行列]
  	for(i=0; i<N_COMPONENT; i++){
		for(k=0; k<N_P; k++){
	  		b[i][k] = 0.0;
	  		for(j=0; j<N_SP; j++){
				b[i][k] += P[j][k] * Val_of_SP[j][i];
				//printf("P[%d][%d]= %le Val=%le\n",j,k,P[j][k],Val_of_SP[j][i] );
	  		}
		}
  	}

  	//係数行列（Pの行列の転置とPを掛けたヤツ。）
  	for(i=0; i<N_P; i++){
		for(j=0; j<N_P; j++){
	  		A[i][j] = 0.0;
	  		for(k=0; k<N_SP; k++){
				A[i][j] += P[k][i] * P[k][j];
	  		}
		}
  	}

  	//逆行列を求める
  	for(i=0; i<N_P; i++){
		for(j=0; j<N_P; j++){
	  		A_Inv[i][j] = A[i][j];
		}
  	}
  	InverseMatrix_4D(A_Inv, ip);

  	//入力した行列と逆行列との積が単位行列になるかチェック
  	for(i=0; i<N_P; i++){
		for(j=0; j<N_P; j++){
	  		check[i][j] = 0.0;
	  		for(k=0; k<N_P; k++){
				check[i][j] += A[i][k] * A_Inv[k][j];
	  		}
		}
  	}
  	for(i=0;i<N_P;i++){
		for(j=0;j<N_P;j++){
			if(i==j)	{
				if(fabs(check[i][j]-1)>0.0000000001) printf("ERROR Failed to create InverseMatrix [%d][%d]!=1  (%f) \n",i,j,check[i][j] );
			}
			else{
				if(fabs(check[i][j])>0.0000000001) printf("ERROR Failed to create InverseMatrix [%d][%d]!=0  (%f)\n",i,j,check[i][j] );
			}
		}	
	}


	//bの左から係数行列[A]の逆行列をかける。この計算結果が局所最小二乗の近似式の係数
  	for(i=0; i<N_COMPONENT; i++){
		for(j=0; j<N_P; j++){
	  		a[i][j] = 0.0;
	  		for(k=0; k<N_P; k++){
				a[i][j] += A_Inv[j][k] * b[i][k];
	  		}	
	  		//printf("iNode=%d a[%d][%d] = %e  \n",ip, i, j, a[i][j]);
		}
  	}

  	//スムージング解用の基底P ※最初にサンプリング点と基準点の差をとって計算しているため、{x,y,z}={0,0,0}を近似式に代入する。
  	P[0][0] = 1.0;
  	for(k=1; k<N_P; k++){
		P[0][k] = 0.0;
  	}

	//スムージング解
  	for(i=0; i<N_COMPONENT; i++){
		if(data==1)JW_on_Node[ip][i] = 0.0;
		else if(data==2)Pai_on_Node[ip][i] = 0.0;
		else if(data==3)Disp_grad_1_on_Node[ip][i] = 0.0;
		else if(data==4)Disp_grad_2_on_Node[ip][i] = 0.0;
		else if(data==5)Disp_grad_3_on_Node[ip][i] = 0.0;
		else if(data==6)W_1_on_Node[ip][i] = 0.0;
		else if(data==7)W_2_on_Node[ip][i] = 0.0;
		else if(data==8)W_3_on_Node[ip][i] = 0.0;
		for(l=0; l<N_P; l++){
	  		if(data==1)JW_on_Node[ip][i] += P[0][l] * a[i][l];
	  		else if(data==2)Pai_on_Node[ip][i] += P[0][l] * a[i][l];
	  		else if(data==3)Disp_grad_1_on_Node[ip][i] += P[0][l] * a[i][l];
	 	 	else if(data==4)Disp_grad_2_on_Node[ip][i] += P[0][l] * a[i][l];
	  		else if(data==5)Disp_grad_3_on_Node[ip][i] += P[0][l] * a[i][l];
	  		else if(data==6)W_1_on_Node[ip][i] += P[0][l] * a[i][l];
	  		else if(data==7)W_2_on_Node[ip][i] += P[0][l] * a[i][l];
	  		else if(data==8)W_3_on_Node[ip][i] += P[0][l] * a[i][l];
		}
  	}
}
//最小二乗近似(二次多項式近似)(P(応力等の数値) = a + bx + cy + dz + exy + fyz + gxz + hxx + iyy + jzz)で近似している。
int LSM_second_basis(int ip, int N_SP, int data, int N_COMPONENT)
{
  	//変数の設定
  	double a[N_COMPONENT][N_P_2]; //近似多項式の係数(これを求める)
  	double P[N_SP][N_P_2]; //基底
  	double A[N_P_2][N_P_2]; //近似多項式をa[i]についての偏微分してできる連立一次方程式におけるaの係数行列
  	double A_Inv[N_P_2][N_P_2]; //行列Aの逆行列
  	double b[N_COMPONENT][N_P_2]; //近似多項式をa[i]についての偏微分してできる連立一次方程式における右辺の項
  	double check[N_P_2][N_P_2];

  	//配列の初期化
  	int i, j, k, l;
  	for(j=0; j<N_SP; j++){
		for(k=0; k<N_P_2; k++){
	  		P[j][k] = 0.0;
		}
  	}
  	for(i=0; i<N_COMPONENT; i++){
		for(k=0; k<N_P_2; k++){
	  		b[i][k] = 0.0;
		}
  	}

  	//基底ベクトルP(二次近似_3次元){1,x,y,z,xy,yz,xz,xx,yy,zz}がサンプリング点数だけある形の行列。
  	for(j=0; j<N_SP; j++){
		P[j][0] = 1.0;
		P[j][1] = Coord_of_SP[j][0] - Coord_of_PatchCenter[ip][0];
		P[j][2] = Coord_of_SP[j][1] - Coord_of_PatchCenter[ip][1];
		P[j][3] = Coord_of_SP[j][2] - Coord_of_PatchCenter[ip][2];
		P[j][4] = (Coord_of_SP[j][0] - Coord_of_PatchCenter[ip][0]) * (Coord_of_SP[j][1] - Coord_of_PatchCenter[ip][1]);
		P[j][5] = (Coord_of_SP[j][1] - Coord_of_PatchCenter[ip][1]) * (Coord_of_SP[j][2] - Coord_of_PatchCenter[ip][2]);
		P[j][6] = (Coord_of_SP[j][2] - Coord_of_PatchCenter[ip][2]) * (Coord_of_SP[j][0] - Coord_of_PatchCenter[ip][0]);
		P[j][7] = (Coord_of_SP[j][0] - Coord_of_PatchCenter[ip][0]) * (Coord_of_SP[j][0] - Coord_of_PatchCenter[ip][0]);
		P[j][8] = (Coord_of_SP[j][1] - Coord_of_PatchCenter[ip][1]) * (Coord_of_SP[j][1] - Coord_of_PatchCenter[ip][1]);
		P[j][9] = (Coord_of_SP[j][2] - Coord_of_PatchCenter[ip][2]) * (Coord_of_SP[j][2] - Coord_of_PatchCenter[ip][2]);
	}

  	//[Pの行列][ガウス点の持つ値の行列]
  	for(i=0; i<N_COMPONENT; i++){
		for(k=0; k<N_P_2; k++){
	  		b[i][k] = 0.0;
	  		for(j=0; j<N_SP; j++){
				b[i][k] += P[j][k] * Val_of_SP[j][i];
	  		}
		}
  	}
  	
  	//係数行列（Pの行列の転置とPを掛けたヤツ。）
  	for(i=0; i<N_P_2; i++){
		for(j=0; j<N_P_2; j++){
	 		A[i][j] = 0.0;
	  		for(k=0; k<N_SP; k++){
				A[i][j] += P[k][i] * P[k][j];
	  		}
		}
  	}

  	//逆行列を求める
  	for(i=0; i<N_P_2; i++){ //行列Aを保持するための措置
		for(j=0; j<N_P_2; j++){
	  		A_Inv[i][j] = A[i][j];
		}
  	}

  	int check_determinant;
	check_determinant=Inverse(A_Inv);

  	////入力した行列と逆行列との積が単位行列になるかチェック
  	for(i=0; i<N_P_2; i++){
		for(j=0; j<N_P_2; j++){
	  		check[i][j] = 0.0;
	  		for(k=0; k<N_P_2; k++){
				check[i][j] += A[i][k] * A_Inv[k][j];
	  		}
		}
  	}
  	for(i=0;i<N_P;i++){
		for(j=0;j<N_P;j++){
			if(i==j)	{
				if(fabs(check[i][j]-1)>0.0000000001) printf("ERROR Failed to create InverseMatrix [%d][%d]!=1  (%f) \n",i,j,check[i][j] );
			}
			else{
				if(fabs(check[i][j])>0.0000000001) printf("ERROR Failed to create InverseMatrix [%d][%d]!=0  (%f)\n",i,j,check[i][j] );
			}
		}	
	}
	//bの左から係数行列[A]の逆行列をかける。この計算結果が局所最小二乗の近似式の係数
	for(i=0; i<N_COMPONENT; i++){
	  	for(j=0; j<N_P_2; j++){
			a[i][j] = 0.0;
			for(k=0; k<N_P_2; k++){
	  			a[i][j] += A_Inv[j][k] * b[i][k];
			}
	  	}
	}

  	//スムージング解用の基底P ※最初にサンプリング点と基準点の差をとって計算しているため、{x,y,z}={0,0,0}を近似式に代入する。
	P[0][0] = 1.0;
	for(k=1; k<N_P_2; k++){
	  	P[0][k] = 0.0;
	}

	//スムージング解
	for(i=0; i<N_COMPONENT; i++){
	  	if(data==1)JW_on_Node[ip][i] = 0.0;
	  	else if(data==2)Pai_on_Node[ip][i] = 0.0;
	  	else if(data==3)Disp_grad_1_on_Node[ip][i] = 0.0;
	  	else if(data==4)Disp_grad_2_on_Node[ip][i] = 0.0;
	  	else if(data==5)Disp_grad_3_on_Node[ip][i] = 0.0;
	  	else if(data==6)W_1_on_Node[ip][i] = 0.0;
	  	else if(data==7)W_2_on_Node[ip][i] = 0.0;
	  	else if(data==8)W_3_on_Node[ip][i] = 0.0;
	  	for(l=0; l<N_P_2; l++){
			if(data==1)JW_on_Node[ip][i] += P[0][l] * a[i][l];
			else if(data==2)Pai_on_Node[ip][i] += P[0][l] * a[i][l];
			else if(data==3)Disp_grad_1_on_Node[ip][i] += P[0][l] * a[i][l];
			else if(data==4)Disp_grad_2_on_Node[ip][i] += P[0][l] * a[i][l];
			else if(data==5)Disp_grad_3_on_Node[ip][i] += P[0][l] * a[i][l];
			else if(data==6)W_1_on_Node[ip][i] += P[0][l] * a[i][l];
			else if(data==7)W_2_on_Node[ip][i] += P[0][l] * a[i][l];
			else if(data==8)W_3_on_Node[ip][i] += P[0][l] * a[i][l];
	  	}	
	}

	return check_determinant;
}
/*****************************************************************************
   積分点のデータを節点に移す。（局所最小二乗(線形多項式、二次多項式)or平均）
*******************************************************************************/
//基準点になる節点を持つ要素を探索、探索後にその要素の積分点から最小二乗近似（上記LSM)をまわす。
int SP_all_Gausspoint_in_eachelement_linear(int data, int N_COMPONENT, int ip)
{
  	int ie, ele, ig, id, in, ic, i, j, k, n, Crack_Front=0;
 	int N_SP, SP;
 	int isp,count,det_check;

	N_SP = 0;
	for(ele=0; ele<NElements_DI; ele++)
	{
		ie = List_Element_DI[ele];
	  	for(i=0; i<N_Vertex_Nodes; i++)//頂点節点のみでLSMを行う。
	  	{
			if(ElementNodeId_s[ie][i] == ip)
			{
		  		for(ig=0; ig<N_IntegrationPoint; ig++)
		  		{
					for(id=0; id<DIMENSION; id++){
		  				Coord_of_SP[N_SP][id] = Gxi_coord[ie][ig][id];//サンプル点の座標
					}
					for(j=0; j<DIMENSION; j++){
		  				for(k=0; k<DIMENSION; k++){
							if(data==1)Val_of_SP[N_SP][0] = JW[ie][ig];
							else if(data==2)Val_of_SP[N_SP][j*3+k] = Pai[ie][ig][j][k];
							else if(data==3)Val_of_SP[N_SP][j] = Disp_grad_1[ie][ig][j];
							else if(data==4)Val_of_SP[N_SP][j] = Disp_grad_2[ie][ig][j];
							else if(data==5)Val_of_SP[N_SP][j] = Disp_grad_3[ie][ig][j];
							else if(data==6)Val_of_SP[N_SP][j] = W_1[ie][ig][j];
							else if(data==7)Val_of_SP[N_SP][j] = W_2[ie][ig][j];
							else if(data==8)Val_of_SP[N_SP][j] = W_3[ie][ig][j];
		  				}
					}
					N_SP++;
	  			}
			}
		}
	}

	for(id=0; id<DIMENSION; id++){
	  	Coord_of_PatchCenter[ip][id] = Node_Coordinate[ip][id];//パッチ中心点の座標
	}
	
  	if(N_SP!=0){
		LSM(ip, N_SP, data, N_COMPONENT);
  	}

	if( isnan(JW_on_Node[ip][0]) != 0 ){
	  	//printf("PCID %d %e %e %e\n", ip, Coord_of_PatchCenter[ip][0], Coord_of_PatchCenter[ip][1], Coord_of_PatchCenter[ip][2]);
	  	for(i=0; i<1; i++){
			//printf(" %e", JW_on_Node[ip][i]);
	  	}
	  	//printf("\n");
	  	for(isp=0; isp<N_SP; isp++){
			//printf("SPID[%d] %e %e %e %e\n", isp, Coord_of_SP[isp][0], Coord_of_SP[isp][1], Coord_of_SP[isp][2], Val_of_SP[isp][0]);
	  	}
	  	//printf("\n");
	}
}

//基準点になる節点を持つ要素を探索、探索後にその要素の積分点から最小二乗近似（上記LSM_second_Basis)をまわす。未整理
int SP_all_Gausspoint_in_eachelement_quad(int data, int N_COMPONENT, int ip)
{
  	int ie, ig, id, in, ic, i, j, k, n, isp, Crack_Front=0;
  	int N_SP, SP;
  	double Gxi_Coord_Ave;

  	for(n=0; n<Total_Crack_Front_Node;n++){
		if(ip==Crack_Front_Node[n]){
	  		Crack_Front=1;
	  		break;
		}
  	}

  	for(isp=0; isp<10000; isp++){
		SP_Element_ID[isp]=0;
  	}

  	N_SP = 0;
  	for(ie=0; ie<NElements; ie++){
		Gxi_Coord_Ave=(Gxi_coord[ie][0][2]+Gxi_coord[ie][1][2]+Gxi_coord[ie][2][2]+Gxi_coord[ie][3][2])/4;
		for(i=0; i<Number_of_Nodes_on_Elements; i++){
	  		if(ElementNodeId_s[ie][i] == ip){
				//printf("%d\n", ie);
				//if(Crack_Front==0 || Gxi_Coord_Ave>5){//き裂前縁でない、もしくはき裂前縁で要素内の重心のz座標が5以上の場合
				for(ig=0; ig<N_IntegrationPoint; ig++){
	  				SP_Element_ID[N_SP]=ie*N_IntegrationPoint+ig;
	  				//printf("SP_Element_ID[i]=%d\tie=%d\n", SP_Element_ID[N_SP], ie);
	  				for(id=0; id<DIMENSION; id++){
						Coord_of_SP[N_SP][id] = Gxi_coord[ie][ig][id];//サンプル点の座標
	  				}//for(id)
	  				for(j=0; j<DIMENSION; j++){
						for(k=0; k<DIMENSION; k++){
		  					if(data==1)Val_of_SP[N_SP][0] = JW[ie][ig];
		 	 				else if(data==2)Val_of_SP[N_SP][j*3+k] = Pai[ie][ig][j][k];
		  					else if(data==3)Val_of_SP[N_SP][j] = Disp_grad_1[ie][ig][j];
		  					else if(data==4)Val_of_SP[N_SP][j] = Disp_grad_2[ie][ig][j];
		  					else if(data==5)Val_of_SP[N_SP][j] = Disp_grad_3[ie][ig][j];
		  					else if(data==6)Val_of_SP[N_SP][j] = W_1[ie][ig][j];
		  					else if(data==7)Val_of_SP[N_SP][j] = W_2[ie][ig][j];
		  					else if(data==8)Val_of_SP[N_SP][j] = W_3[ie][ig][j];
						}
	  				}
	  				N_SP++;
				}
	  		}
		}
  	}
  	for(id=0; id<3; id++){
		Coord_of_PatchCenter[ip][id] = Node_Coordinate[ip][id];//パッチ中心点の座標
  	}

  	int count;
 	int det_check;
  	int max_N_SP=10;

  	if(N_SP<max_N_SP){
		SP=N_SP;
		N_SP=max_N_SP;
		N_SP=Get_SP_in_circle_distance(data, N_COMPONENT, ip, N_SP, SP, 5);
  	}

  	det_check=0;
  	det_check=LSM_second_basis(ip, N_SP, data, N_COMPONENT);

  	while(det_check==0){
		if(det_check==0){
	  		SP=N_SP;
	  		if(N_SP<100)N_SP+=1;
	  		else if(100<=N_SP && N_SP<200)N_SP+=2;
	  		else N_SP+=10;
	  		//printf("N_SP=%d, SP=%d\n", N_SP, SP);
	  		N_SP=Get_SP_in_circle_distance(data, N_COMPONENT, ip, N_SP, SP, 5);
		}
		det_check=LSM_second_basis(ip, N_SP, data, N_COMPONENT);
		//if( det_check==1){
	  		//printf("det_check=%d\n", det_check);
	  		//printf("node=%d N_SP=%d data=%d  %e\n", ip, N_SP, data, JW_on_Node[ip][0]);
	 		//printf("%lf %lf %lf\n", Coord_of_PatchCenter[ip][0], Coord_of_PatchCenter[ip][1], Coord_of_PatchCenter[ip][2]);
	  		//for(isp=0; isp<N_SP; isp++){
				//printf("SPID[%d] %e %e %e %e\n", isp, Coord_of_SP[isp][0], Coord_of_SP[isp][1], Coord_of_SP[isp][2], Val_of_SP[isp][0]);
				//printf("%e %e %e\n", Coord_of_SP[isp][0], Coord_of_SP[isp][1], Coord_of_SP[isp][2]);
				//printf("%e\n", Val_of_SP[isp][0]);
				//printf("%le %le %le \n", Coord_of_SP[isp][0]-Coord_of_PatchCenter[ip][0], Coord_of_SP[isp][1]-Coord_of_PatchCenter[ip][1], Coord_of_SP[isp][2]-Coord_of_PatchCenter[ip][2]);
	  		//}
		//}
  	}

	double check_data=0.0;
	if(data==1)check_data=JW_on_Node[ip][0];
	if( isnan(check_data) != 0 || check_data<0.0){
	  	//printf("node=%d N_SP=%d data=%d  %e\n", ip, N_SP, data, JW_on_Node[ip][0]);
	  	//printf("%lf %lf %lf\n", Coord_of_PatchCenter[ip][0], Coord_of_PatchCenter[ip][1], Coord_of_PatchCenter[ip][2]);
	  	//for(isp=0; isp<N_SP; isp++){
			//printf("SPID[%d] %e %e %e %e\n", isp, Coord_of_SP[isp][0], Coord_of_SP[isp][1], Coord_of_SP[isp][2], Val_of_SP[isp][0]);
			//printf("%e %e %e\n", Coord_of_SP[isp][0], Coord_of_SP[isp][1], Coord_of_SP[isp][2]);
			//printf("%e\n", Val_of_SP[isp][0]);
			//printf("%le %le %le \n", Coord_of_SP[isp][0]-Coord_of_PatchCenter[ip][0], Coord_of_SP[isp][1]-Coord_of_PatchCenter[ip][1], Coord_of_SP[isp][2]-Coord_of_PatchCenter[ip][2]);
	  	//}
	}

  	for(isp=0; isp<N_SP; isp++){
		SP_Element_ID[isp]=0;
  	}
}

//未整理、上の〜〜〜Quadで使ってるっぽい。
int Get_SP_in_circle_distance(int data, int N_COMPONENT, int ip, int N_SP, int SP, double box)
{
  	int ie, ig, id, in, ic, i,j,k,isp, dummy_ID, dummy_in, SPID[N_SP], count, flag_num=0, dummy_IE;
  	double r, r2, d, c1=0.5, c2=1, distance[N_SP], dummy;
  	//printf("box=%le\nthreshold=%le\n", box, threshold);
  	//for(isp=0; isp<N_SP; isp++){
  	//  SP_Element_ID[isp]=0;
  	//}

  	for(isp=SP; isp<N_SP; isp++){
		distance[isp] = 1000000.0;
		SPID[isp]=-1;
  	}

  	if(Crack_flag[ip]==1)flag_num=2;
  	else if(Crack_flag[ip]==2)flag_num=1;

	//printf("flag_num=%d\n", flag_num);
  	for(i=0; i<SP; i++){
		//printf("SP_Element_ID[%d]=%d\n", i, SP_Element_ID[i]);
  	}
  	for(ie=0; ie<NElements; ie++){
		for(ig=0; ig<N_IntegrationPoint; ig++){
	  		if(   Node_Coordinate[ip][0]-box <= Gxi_coord[ie][ig][0] && Gxi_coord[ie][ig][0] <= Node_Coordinate[ip][0]+box
		 		&& Node_Coordinate[ip][1]-box <= Gxi_coord[ie][ig][1] && Gxi_coord[ie][ig][1] <= Node_Coordinate[ip][1]+box
		 		&& Node_Coordinate[ip][2]-box <= Gxi_coord[ie][ig][2] && Gxi_coord[ie][ig][2] <= Node_Coordinate[ip][2]+box ){
				
				dummy_in=ie*N_IntegrationPoint+ig;
				count=0;
				for(i=0; i<SP; i++){
	  				if(SP_Element_ID[i]==ie*N_IntegrationPoint+ig){
						//printf("SP_Element_ID[i]=%d\tie*4+ig=%d\n", SP_Element_ID[i], ie*4+ig);
						count++;
	  				}
				}
				//printf("ip=%d\tcount=%d\n", ip, count);
				if(count==0 && Gauss_flag[ie]!=flag_num){
	  				r2 = 0.0;
	  				for(id=0; id<DIMENSION; id++){
						r2 += pow(Gxi_coord[ie][ig][id] - Node_Coordinate[ip][id], 2);
	  				}
	  				r = sqrt(r2);
	  				for(isp=SP; isp<N_SP; isp++){
						if (r <= distance[isp]){
		  					dummy = distance[isp];
		  					distance[isp] = r;
		  					r = dummy;	
		  					dummy_ID = SPID[isp];
		  					SPID[isp] = dummy_in;
		  					dummy_in = dummy_ID;
		  					dummy_IE = ie;
						}//if(distance)
	  				}//for(isp)
				}//if(count)
	  		}//if(box)
		}//for(ig)
  	}//for(e)

 	//printf("Gauss_flag[%d]=%d\n", dummy_IE, Gauss_flag[dummy_IE]);
  	//for(id=0; id<3; id++){
  	//  Coord_of_PatchCenter[ip][id] = Node_Coordinate[ip][id];//パッチ中心点の座標
  	//}
  	for(isp=SP; isp<N_SP; isp++){
		//printf("SPID[%d]=%d\n", isp, SPID[isp]);
		//printf("distance[%d]=%lf\n", isp, distance[isp]);
  	}
	//SPID[0]=-1;
  	for(isp=SP; isp<N_SP; isp++){
		if(SPID[isp]==-1){
	  		//printf("isp=%d\n", isp);
	  		//N_SP=isp;
	  		//break;
			printf("ERROR ip=%d isp=%d is not enough sample points\n", ip, isp);
			exit (1);
		}
  	}
  	//printf("N_SP=%d\n", N_SP);

  	for(isp=SP; isp<N_SP; isp++){
		SP_Element_ID[isp]=SPID[isp];
		for(j=0; j<DIMENSION; j++){
	  		for(k=0; k<DIMENSION; k++){
				for(id=0; id<DIMENSION; id++){
	  				Coord_of_SP[isp][id] = Gxi_coord[SPID[isp]/4][SPID[isp]%4][id];//サンプル点での物理量
		  			//printf("element=%d GaussPoint=%d\n", SPID[isp]/4, SPID[isp]%4);
				}
				if(data==1)Val_of_SP[isp][0] = JW[SPID[isp]/4][SPID[isp]%4];
				//if(data==1)Val_of_SP[isp][0] = Gxi_coord[SPID[isp]/4][SPID[isp]%4][0]*Gxi_coord[SPID[isp]/4][SPID[isp]%4][0];
				else if(data==2)Val_of_SP[isp][j*3+k] = Pai[SPID[isp]/4][SPID[isp]%4][j][k];
				//if(data==2)Val_of_SP[isp][j*3+k] = Gxi_coord[SPID[isp]/4][SPID[isp]%4][0];
				else if(data==3)Val_of_SP[isp][j] = Disp_grad_1[SPID[isp]/4][SPID[isp]%4][j];
				//if(data==3)Val_of_SP[isp][j] = Gxi_coord[SPID[isp]/4][SPID[isp]%4][0];
				else if(data==4)Val_of_SP[isp][j] = Disp_grad_2[SPID[isp]/4][SPID[isp]%4][j];
				//if(data==4)Val_of_SP[isp][j] = Gxi_coord[SPID[isp]/4][SPID[isp]%4][0];
				else if(data==5)Val_of_SP[isp][j] = Disp_grad_3[SPID[isp]/4][SPID[isp]%4][j];
				//if(data==5)Val_of_SP[isp][j] = Gxi_coord[SPID[isp]/4][SPID[isp]%4][0];
				else if(data==6)Val_of_SP[isp][j] = W_1[SPID[isp]/4][SPID[isp]%4][j];
				//if(data==6)Val_of_SP[isp][j] = Gxi_coord[SPID[isp]/4][SPID[isp]%4][0];
				else if(data==7)Val_of_SP[isp][j] = W_2[SPID[isp]/4][SPID[isp]%4][j];
				//if(data==6)Val_of_SP[isp][j] = Gxi_coord[SPID[isp]/4][SPID[isp]%4][0];
				else if(data==8)Val_of_SP[isp][j] = W_3[SPID[isp]/4][SPID[isp]%4][j];
				//if(data==6)Val_of_SP[isp][j] = Gxi_coord[SPID[isp]/4][SPID[isp]%4][0];
	  		}//for(k)
	  		//Val_of_SP[N_SP][0] = Gxi_coord[ie][ig][0]*10;
		}//for(j)
  	}
  	//printf("N_SP=%d\n", N_SP);
  	return N_SP;
}

//基準点にする節点の周りのガウス点から普通に平均とる。
int SP_Average_one_Gausspoint_in_eachelement(int data, int N_COMPONENT, int ip)
{
  	int ie, ig, id, in, ic, i,j,k;
  	int N_SP;
 	int l_min;
  	N_SP = 0;
  	int isp;
  	int count;

	//初期化？
   	for(j=0; j<DIMENSION; j++){
	 	for(k=0; k<DIMENSION; k++){
	   		JW_on_Node[ip][0]=0.0;
	   		Pai_on_Node[ip][j*3+k]=0.0;
	   		Disp_grad_1_on_Node[ip][j]=0.0;
	   		Disp_grad_2_on_Node[ip][j]=0.0;
	   		Disp_grad_3_on_Node[ip][j]=0.0;
	   		W_1_on_Node[ip][j]=0.0;
	   		W_2_on_Node[ip][j]=0.0;
	   		W_3_on_Node[ip][j]=0.0;
	 	}
   	}

  	for(ie=0; ie<NElements; ie++){
		for(i=0; i<N_Vertex_Nodes; i++){
	  		if(ElementNodeId_s[ie][i] == ip){
				l_min=0;
				for(id=0; id<DIMENSION; id++){
	  				Coord_of_SP[N_SP][id] = Gxi_coord[ie][0][id];//サンプル点の座標
				}
				for(ig=1; ig<N_IntegrationPoint; ig++){
	  				if((pow((Gxi_coord[ie][ig][0]-Node_Coordinate[ip][0]),2)
		  				+pow((Gxi_coord[ie][ig][1]-Node_Coordinate[ip][1]),2)
		  				+pow((Gxi_coord[ie][ig][2]-Node_Coordinate[ip][2]),2))
		 			< (pow((Coord_of_SP[N_SP][0]-Node_Coordinate[ip][0]),2)
						+pow((Coord_of_SP[N_SP][1]-Node_Coordinate[ip][1]),2)
						+pow((Coord_of_SP[N_SP][2]-Node_Coordinate[ip][2]),2))){
						l_min=ig;
						for(id=0; id<DIMENSION; id++){
		  					Coord_of_SP[N_SP][id] = Gxi_coord[ie][ig][id];//サンプル点の座標(更新後)
						}
	  				}
				}

				if(data==1)Val_of_SP[N_SP][0] = JW[ie][l_min];
				else{
					for(j=0; j<DIMENSION; j++){
						if(data==2){
							for(k=0; k<DIMENSION; k++){
								Val_of_SP[N_SP][j*3+k] = Pai[ie][l_min][j][k];
	  						}
						}
						else if(data==3)Val_of_SP[N_SP][j] = Disp_grad_1[ie][l_min][j];
						else if(data==4)Val_of_SP[N_SP][j] = Disp_grad_2[ie][l_min][j];
						else if(data==5)Val_of_SP[N_SP][j] = Disp_grad_3[ie][l_min][j];
						else if(data==6)Val_of_SP[N_SP][j] = W_1[ie][l_min][j];
						else if(data==7)Val_of_SP[N_SP][j] = W_2[ie][l_min][j];
						else if(data==8)Val_of_SP[N_SP][j] = W_3[ie][l_min][j];
					}
				}
				N_SP++;
	  		}
		}
	}
  	
  	for(id=0; id<DIMENSION; id++){
		Coord_of_PatchCenter[ip][id] = Node_Coordinate[ip][id];//パッチ中心点の座標
  	}

 	for(isp=0; isp<N_SP; isp++){
   		if(data==1)JW_on_Node[ip][0] += Val_of_SP[isp][0];
   		else if(data==2){
   			for(ic=0; ic<9; ic++){
	  			Pai_on_Node[ip][ic] += Val_of_SP[isp][ic];
			}
   		}
   		else{
   			for(j=0; j<DIMENSION; j++){
	  			if(data==3)Disp_grad_1_on_Node[ip][j] += Val_of_SP[isp][j];
	  			else if(data==4)Disp_grad_2_on_Node[ip][j] += Val_of_SP[isp][j];
	  			else if(data==5)Disp_grad_3_on_Node[ip][j] += Val_of_SP[isp][j];
	  			else if(data==6)W_1_on_Node[ip][j] += Val_of_SP[isp][j];
	  			else if(data==7)W_2_on_Node[ip][j] += Val_of_SP[isp][j];
	  			else if(data==8)W_3_on_Node[ip][j] += Val_of_SP[isp][j];
			}
   		}

  	}

   	if(data==1)JW_on_Node[ip][0] = JW_on_Node[ip][0]/N_SP;
   	else if(data==2)
		for(ic=0; ic<9; ic++){
	  		Pai_on_Node[ip][ic] = Pai_on_Node[ip][ic]/N_SP;
		}
	else{
		for(j=0; j<DIMENSION; j++){
	  		if(data==3)Disp_grad_1_on_Node[ip][j] = Disp_grad_1_on_Node[ip][j]/N_SP;
	  		else if(data==4)Disp_grad_2_on_Node[ip][j] = Disp_grad_2_on_Node[ip][j]/N_SP;
	  		else if(data==5)Disp_grad_3_on_Node[ip][j] = Disp_grad_3_on_Node[ip][j]/N_SP;
	  		else if(data==6)W_1_on_Node[ip][j] = W_1_on_Node[ip][j]/N_SP;
	  		else if(data==7)W_2_on_Node[ip][j] = W_2_on_Node[ip][j]/N_SP;
	  		else if(data==8)W_3_on_Node[ip][j] = W_3_on_Node[ip][j]/N_SP;
		}
	}
}

/*****************************************************************************
   			二次要素時に中間節点に平均化したデータを突っ込む
*******************************************************************************/
//もともとIn house形式のコネクティビティだったものを.msh形式のコネクティビティに変更。T10専用、頂点節点のデータの平均値を中間節点に格納する
void data_to_intermediate_node_T10()
{
  	int ip, i, ie, j, k;

  	for(ie=0; ie<NElements; ie++){
		//for(i=0; i<NO_NODES_ON_ELEMENT; i++){
		JW_on_Node[ElementNodeId_s[ie][4]][0] = 0.5*(JW_on_Node[ElementNodeId_s[ie][0]][0] + JW_on_Node[ElementNodeId_s[ie][1]][0]);
		JW_on_Node[ElementNodeId_s[ie][5]][0] = 0.5*(JW_on_Node[ElementNodeId_s[ie][0]][0] + JW_on_Node[ElementNodeId_s[ie][2]][0]);
		JW_on_Node[ElementNodeId_s[ie][6]][0] = 0.5*(JW_on_Node[ElementNodeId_s[ie][0]][0] + JW_on_Node[ElementNodeId_s[ie][3]][0]);
		JW_on_Node[ElementNodeId_s[ie][7]][0] = 0.5*(JW_on_Node[ElementNodeId_s[ie][1]][0] + JW_on_Node[ElementNodeId_s[ie][2]][0]);
		JW_on_Node[ElementNodeId_s[ie][8]][0] = 0.5*(JW_on_Node[ElementNodeId_s[ie][2]][0] + JW_on_Node[ElementNodeId_s[ie][3]][0]);
		JW_on_Node[ElementNodeId_s[ie][9]][0] = 0.5*(JW_on_Node[ElementNodeId_s[ie][3]][0] + JW_on_Node[ElementNodeId_s[ie][1]][0]);

		//printf("%le %le %le %le \n", JW_on_Node[ElementNodeId_s[ie][0]][0], JW_on_Node[ElementNodeId_s[ie][1]][0], JW_on_Node[ElementNodeId_s[ie][2]][0], JW_on_Node[ElementNodeId_s[ie][3]][0]);
		//printf("%le %le %le \n", JW_on_Node[ElementNodeId_s[ie][4]][0], JW_on_Node[ElementNodeId_s[ie][0]][0], JW_on_Node[ElementNodeId_s[ie][1]][0]);

	  	for(j=0; j<3; j++){
			for(k=0; k<3; k++){
				Pai_on_Node[ElementNodeId_s[ie][4]][j*3+k] = 0.5*(Pai_on_Node[ElementNodeId_s[ie][0]][j*3+k] + Pai_on_Node[ElementNodeId_s[ie][1]][j*3+k]);
				Pai_on_Node[ElementNodeId_s[ie][5]][j*3+k] = 0.5*(Pai_on_Node[ElementNodeId_s[ie][0]][j*3+k] + Pai_on_Node[ElementNodeId_s[ie][2]][j*3+k]);
				Pai_on_Node[ElementNodeId_s[ie][6]][j*3+k] = 0.5*(Pai_on_Node[ElementNodeId_s[ie][0]][j*3+k] + Pai_on_Node[ElementNodeId_s[ie][3]][j*3+k]);
				Pai_on_Node[ElementNodeId_s[ie][7]][j*3+k] = 0.5*(Pai_on_Node[ElementNodeId_s[ie][1]][j*3+k] + Pai_on_Node[ElementNodeId_s[ie][2]][j*3+k]);
				Pai_on_Node[ElementNodeId_s[ie][8]][j*3+k] = 0.5*(Pai_on_Node[ElementNodeId_s[ie][2]][j*3+k] + Pai_on_Node[ElementNodeId_s[ie][3]][j*3+k]);
				Pai_on_Node[ElementNodeId_s[ie][9]][j*3+k] = 0.5*(Pai_on_Node[ElementNodeId_s[ie][3]][j*3+k] + Pai_on_Node[ElementNodeId_s[ie][1]][j*3+k]);
		 	}

			Disp_grad_1_on_Node[ElementNodeId_s[ie][4]][j] = 0.5*(Disp_grad_1_on_Node[ElementNodeId_s[ie][0]][j] + Disp_grad_1_on_Node[ElementNodeId_s[ie][1]][j]);
			Disp_grad_1_on_Node[ElementNodeId_s[ie][5]][j] = 0.5*(Disp_grad_1_on_Node[ElementNodeId_s[ie][0]][j] + Disp_grad_1_on_Node[ElementNodeId_s[ie][2]][j]);
			Disp_grad_1_on_Node[ElementNodeId_s[ie][6]][j] = 0.5*(Disp_grad_1_on_Node[ElementNodeId_s[ie][0]][j] + Disp_grad_1_on_Node[ElementNodeId_s[ie][3]][j]);
			Disp_grad_1_on_Node[ElementNodeId_s[ie][7]][j] = 0.5*(Disp_grad_1_on_Node[ElementNodeId_s[ie][1]][j] + Disp_grad_1_on_Node[ElementNodeId_s[ie][2]][j]);
			Disp_grad_1_on_Node[ElementNodeId_s[ie][8]][j] = 0.5*(Disp_grad_1_on_Node[ElementNodeId_s[ie][2]][j] + Disp_grad_1_on_Node[ElementNodeId_s[ie][3]][j]);
			Disp_grad_1_on_Node[ElementNodeId_s[ie][9]][j] = 0.5*(Disp_grad_1_on_Node[ElementNodeId_s[ie][3]][j] + Disp_grad_1_on_Node[ElementNodeId_s[ie][1]][j]);

			Disp_grad_2_on_Node[ElementNodeId_s[ie][4]][j] = 0.5*(Disp_grad_2_on_Node[ElementNodeId_s[ie][0]][j] + Disp_grad_2_on_Node[ElementNodeId_s[ie][1]][j]);
			Disp_grad_2_on_Node[ElementNodeId_s[ie][5]][j] = 0.5*(Disp_grad_2_on_Node[ElementNodeId_s[ie][0]][j] + Disp_grad_2_on_Node[ElementNodeId_s[ie][2]][j]);
			Disp_grad_2_on_Node[ElementNodeId_s[ie][6]][j] = 0.5*(Disp_grad_2_on_Node[ElementNodeId_s[ie][0]][j] + Disp_grad_2_on_Node[ElementNodeId_s[ie][3]][j]);
			Disp_grad_2_on_Node[ElementNodeId_s[ie][7]][j] = 0.5*(Disp_grad_2_on_Node[ElementNodeId_s[ie][1]][j] + Disp_grad_2_on_Node[ElementNodeId_s[ie][2]][j]);
			Disp_grad_2_on_Node[ElementNodeId_s[ie][8]][j] = 0.5*(Disp_grad_2_on_Node[ElementNodeId_s[ie][2]][j] + Disp_grad_2_on_Node[ElementNodeId_s[ie][3]][j]);
			Disp_grad_2_on_Node[ElementNodeId_s[ie][9]][j] = 0.5*(Disp_grad_2_on_Node[ElementNodeId_s[ie][3]][j] + Disp_grad_2_on_Node[ElementNodeId_s[ie][1]][j]);

			Disp_grad_3_on_Node[ElementNodeId_s[ie][4]][j] = 0.5*(Disp_grad_3_on_Node[ElementNodeId_s[ie][0]][j] + Disp_grad_3_on_Node[ElementNodeId_s[ie][1]][j]);
			Disp_grad_3_on_Node[ElementNodeId_s[ie][5]][j] = 0.5*(Disp_grad_3_on_Node[ElementNodeId_s[ie][0]][j] + Disp_grad_3_on_Node[ElementNodeId_s[ie][2]][j]);
			Disp_grad_3_on_Node[ElementNodeId_s[ie][6]][j] = 0.5*(Disp_grad_3_on_Node[ElementNodeId_s[ie][0]][j] + Disp_grad_3_on_Node[ElementNodeId_s[ie][3]][j]);
			Disp_grad_3_on_Node[ElementNodeId_s[ie][7]][j] = 0.5*(Disp_grad_3_on_Node[ElementNodeId_s[ie][1]][j] + Disp_grad_3_on_Node[ElementNodeId_s[ie][2]][j]);
			Disp_grad_3_on_Node[ElementNodeId_s[ie][8]][j] = 0.5*(Disp_grad_3_on_Node[ElementNodeId_s[ie][2]][j] + Disp_grad_3_on_Node[ElementNodeId_s[ie][3]][j]);
			Disp_grad_3_on_Node[ElementNodeId_s[ie][9]][j] = 0.5*(Disp_grad_3_on_Node[ElementNodeId_s[ie][3]][j] + Disp_grad_3_on_Node[ElementNodeId_s[ie][1]][j]);

			W_1_on_Node[ElementNodeId_s[ie][4]][j] = 0.5*(W_1_on_Node[ElementNodeId_s[ie][0]][j] + W_1_on_Node[ElementNodeId_s[ie][1]][j]);
			W_1_on_Node[ElementNodeId_s[ie][5]][j] = 0.5*(W_1_on_Node[ElementNodeId_s[ie][0]][j] + W_1_on_Node[ElementNodeId_s[ie][2]][j]);
			W_1_on_Node[ElementNodeId_s[ie][6]][j] = 0.5*(W_1_on_Node[ElementNodeId_s[ie][0]][j] + W_1_on_Node[ElementNodeId_s[ie][3]][j]);
			W_1_on_Node[ElementNodeId_s[ie][7]][j] = 0.5*(W_1_on_Node[ElementNodeId_s[ie][1]][j] + W_1_on_Node[ElementNodeId_s[ie][2]][j]);
			W_1_on_Node[ElementNodeId_s[ie][8]][j] = 0.5*(W_1_on_Node[ElementNodeId_s[ie][2]][j] + W_1_on_Node[ElementNodeId_s[ie][3]][j]);
			W_1_on_Node[ElementNodeId_s[ie][9]][j] = 0.5*(W_1_on_Node[ElementNodeId_s[ie][3]][j] + W_1_on_Node[ElementNodeId_s[ie][1]][j]);

			W_2_on_Node[ElementNodeId_s[ie][4]][j] = 0.5*(W_2_on_Node[ElementNodeId_s[ie][0]][j] + W_2_on_Node[ElementNodeId_s[ie][1]][j]);
			W_2_on_Node[ElementNodeId_s[ie][5]][j] = 0.5*(W_2_on_Node[ElementNodeId_s[ie][0]][j] + W_2_on_Node[ElementNodeId_s[ie][2]][j]);
			W_2_on_Node[ElementNodeId_s[ie][6]][j] = 0.5*(W_2_on_Node[ElementNodeId_s[ie][0]][j] + W_2_on_Node[ElementNodeId_s[ie][3]][j]);
			W_2_on_Node[ElementNodeId_s[ie][7]][j] = 0.5*(W_2_on_Node[ElementNodeId_s[ie][1]][j] + W_2_on_Node[ElementNodeId_s[ie][2]][j]);
			W_2_on_Node[ElementNodeId_s[ie][8]][j] = 0.5*(W_2_on_Node[ElementNodeId_s[ie][2]][j] + W_2_on_Node[ElementNodeId_s[ie][3]][j]);
			W_2_on_Node[ElementNodeId_s[ie][9]][j] = 0.5*(W_2_on_Node[ElementNodeId_s[ie][3]][j] + W_2_on_Node[ElementNodeId_s[ie][1]][j]);

			W_3_on_Node[ElementNodeId_s[ie][4]][j] = 0.5*(W_3_on_Node[ElementNodeId_s[ie][0]][j] + W_3_on_Node[ElementNodeId_s[ie][1]][j]);
			W_3_on_Node[ElementNodeId_s[ie][5]][j] = 0.5*(W_3_on_Node[ElementNodeId_s[ie][0]][j] + W_3_on_Node[ElementNodeId_s[ie][2]][j]);
			W_3_on_Node[ElementNodeId_s[ie][6]][j] = 0.5*(W_3_on_Node[ElementNodeId_s[ie][0]][j] + W_3_on_Node[ElementNodeId_s[ie][3]][j]);
			W_3_on_Node[ElementNodeId_s[ie][7]][j] = 0.5*(W_3_on_Node[ElementNodeId_s[ie][1]][j] + W_3_on_Node[ElementNodeId_s[ie][2]][j]);
			W_3_on_Node[ElementNodeId_s[ie][8]][j] = 0.5*(W_3_on_Node[ElementNodeId_s[ie][2]][j] + W_3_on_Node[ElementNodeId_s[ie][3]][j]);
			W_3_on_Node[ElementNodeId_s[ie][9]][j] = 0.5*(W_3_on_Node[ElementNodeId_s[ie][3]][j] + W_3_on_Node[ElementNodeId_s[ie][1]][j]);
	   	}
  	}

  	for(ip=0; ip<NNodes; ip++){
  		if( isnan(JW_on_Node[ip][0]) != 0 ){
	  		//printf("%d %d\n",ip, N_SP);
	  		//printf("PCID %d %e %e %e\n", ip, Coord_of_PatchCenter[ip][0], Coord_of_PatchCenter[ip][1], Coord_of_PatchCenter[ip][2]);
	  		//printf(" %e %e %e", Coord_of_PatchCenter[ip][0], Coord_of_PatchCenter[ip][1], Coord_of_PatchCenter[ip][2]);
	  		for(i=0; i<1; i++){
				//printf(" %e", JW_on_Node[ip][i]);
	  		}
	  		//printf("\n");
		}
  	}
}
//H20専用、頂点節点のデータの平均値を中間節点に格納する
void data_to_intermediate_node_H20()
{
  	int ip, i, ie, j, k;

  	for(ie=0; ie<NElements; ie++){
		JW_on_Node[ElementNodeId_s[ie][8]][0]  = 0.5*(JW_on_Node[ElementNodeId_s[ie][0]][0] + JW_on_Node[ElementNodeId_s[ie][1]][0]);
		JW_on_Node[ElementNodeId_s[ie][9]][0]  = 0.5*(JW_on_Node[ElementNodeId_s[ie][1]][0] + JW_on_Node[ElementNodeId_s[ie][2]][0]);
		JW_on_Node[ElementNodeId_s[ie][10]][0] = 0.5*(JW_on_Node[ElementNodeId_s[ie][2]][0] + JW_on_Node[ElementNodeId_s[ie][3]][0]);
		JW_on_Node[ElementNodeId_s[ie][11]][0] = 0.5*(JW_on_Node[ElementNodeId_s[ie][3]][0] + JW_on_Node[ElementNodeId_s[ie][0]][0]);

		JW_on_Node[ElementNodeId_s[ie][12]][0] = 0.5*(JW_on_Node[ElementNodeId_s[ie][0]][0] + JW_on_Node[ElementNodeId_s[ie][4]][0]);
		JW_on_Node[ElementNodeId_s[ie][13]][0] = 0.5*(JW_on_Node[ElementNodeId_s[ie][1]][0] + JW_on_Node[ElementNodeId_s[ie][5]][0]);
		JW_on_Node[ElementNodeId_s[ie][14]][0] = 0.5*(JW_on_Node[ElementNodeId_s[ie][2]][0] + JW_on_Node[ElementNodeId_s[ie][6]][0]);
		JW_on_Node[ElementNodeId_s[ie][15]][0] = 0.5*(JW_on_Node[ElementNodeId_s[ie][3]][0] + JW_on_Node[ElementNodeId_s[ie][7]][0]);

		JW_on_Node[ElementNodeId_s[ie][16]][0] = 0.5*(JW_on_Node[ElementNodeId_s[ie][4]][0] + JW_on_Node[ElementNodeId_s[ie][5]][0]);
		JW_on_Node[ElementNodeId_s[ie][17]][0] = 0.5*(JW_on_Node[ElementNodeId_s[ie][5]][0] + JW_on_Node[ElementNodeId_s[ie][6]][0]);
		JW_on_Node[ElementNodeId_s[ie][18]][0] = 0.5*(JW_on_Node[ElementNodeId_s[ie][6]][0] + JW_on_Node[ElementNodeId_s[ie][7]][0]);
		JW_on_Node[ElementNodeId_s[ie][19]][0] = 0.5*(JW_on_Node[ElementNodeId_s[ie][7]][0] + JW_on_Node[ElementNodeId_s[ie][4]][0]);

	  	for(j=0; j<3; j++){
			for(k=0; k<3; k++){
				Pai_on_Node[ElementNodeId_s[ie][8]][j*3+k]  = 0.5*(Pai_on_Node[ElementNodeId_s[ie][0]][j*3+k] + Pai_on_Node[ElementNodeId_s[ie][1]][j*3+k]);
				Pai_on_Node[ElementNodeId_s[ie][9]][j*3+k]  = 0.5*(Pai_on_Node[ElementNodeId_s[ie][1]][j*3+k] + Pai_on_Node[ElementNodeId_s[ie][2]][j*3+k]);
				Pai_on_Node[ElementNodeId_s[ie][10]][j*3+k] = 0.5*(Pai_on_Node[ElementNodeId_s[ie][2]][j*3+k] + Pai_on_Node[ElementNodeId_s[ie][3]][j*3+k]);
				Pai_on_Node[ElementNodeId_s[ie][11]][j*3+k] = 0.5*(Pai_on_Node[ElementNodeId_s[ie][3]][j*3+k] + Pai_on_Node[ElementNodeId_s[ie][0]][j*3+k]);
				Pai_on_Node[ElementNodeId_s[ie][12]][j*3+k] = 0.5*(Pai_on_Node[ElementNodeId_s[ie][0]][j*3+k] + Pai_on_Node[ElementNodeId_s[ie][4]][j*3+k]);
				Pai_on_Node[ElementNodeId_s[ie][13]][j*3+k] = 0.5*(Pai_on_Node[ElementNodeId_s[ie][1]][j*3+k] + Pai_on_Node[ElementNodeId_s[ie][5]][j*3+k]);
				Pai_on_Node[ElementNodeId_s[ie][14]][j*3+k] = 0.5*(Pai_on_Node[ElementNodeId_s[ie][2]][j*3+k] + Pai_on_Node[ElementNodeId_s[ie][6]][j*3+k]);
				Pai_on_Node[ElementNodeId_s[ie][15]][j*3+k] = 0.5*(Pai_on_Node[ElementNodeId_s[ie][3]][j*3+k] + Pai_on_Node[ElementNodeId_s[ie][7]][j*3+k]);
				Pai_on_Node[ElementNodeId_s[ie][16]][j*3+k] = 0.5*(Pai_on_Node[ElementNodeId_s[ie][4]][j*3+k] + Pai_on_Node[ElementNodeId_s[ie][5]][j*3+k]);
				Pai_on_Node[ElementNodeId_s[ie][17]][j*3+k] = 0.5*(Pai_on_Node[ElementNodeId_s[ie][5]][j*3+k] + Pai_on_Node[ElementNodeId_s[ie][6]][j*3+k]);
				Pai_on_Node[ElementNodeId_s[ie][18]][j*3+k] = 0.5*(Pai_on_Node[ElementNodeId_s[ie][6]][j*3+k] + Pai_on_Node[ElementNodeId_s[ie][7]][j*3+k]);
				Pai_on_Node[ElementNodeId_s[ie][19]][j*3+k] = 0.5*(Pai_on_Node[ElementNodeId_s[ie][7]][j*3+k] + Pai_on_Node[ElementNodeId_s[ie][4]][j*3+k]);
		 	}

			Disp_grad_1_on_Node[ElementNodeId_s[ie][8]][j]  = 0.5*(Disp_grad_1_on_Node[ElementNodeId_s[ie][0]][j] + Disp_grad_1_on_Node[ElementNodeId_s[ie][1]][j]);
			Disp_grad_1_on_Node[ElementNodeId_s[ie][9]][j]  = 0.5*(Disp_grad_1_on_Node[ElementNodeId_s[ie][1]][j] + Disp_grad_1_on_Node[ElementNodeId_s[ie][2]][j]);
			Disp_grad_1_on_Node[ElementNodeId_s[ie][10]][j] = 0.5*(Disp_grad_1_on_Node[ElementNodeId_s[ie][2]][j] + Disp_grad_1_on_Node[ElementNodeId_s[ie][3]][j]);
			Disp_grad_1_on_Node[ElementNodeId_s[ie][11]][j] = 0.5*(Disp_grad_1_on_Node[ElementNodeId_s[ie][3]][j] + Disp_grad_1_on_Node[ElementNodeId_s[ie][0]][j]);
			Disp_grad_1_on_Node[ElementNodeId_s[ie][12]][j] = 0.5*(Disp_grad_1_on_Node[ElementNodeId_s[ie][0]][j] + Disp_grad_1_on_Node[ElementNodeId_s[ie][4]][j]);
			Disp_grad_1_on_Node[ElementNodeId_s[ie][13]][j] = 0.5*(Disp_grad_1_on_Node[ElementNodeId_s[ie][1]][j] + Disp_grad_1_on_Node[ElementNodeId_s[ie][5]][j]);
			Disp_grad_1_on_Node[ElementNodeId_s[ie][14]][j] = 0.5*(Disp_grad_1_on_Node[ElementNodeId_s[ie][2]][j] + Disp_grad_1_on_Node[ElementNodeId_s[ie][6]][j]);
			Disp_grad_1_on_Node[ElementNodeId_s[ie][15]][j] = 0.5*(Disp_grad_1_on_Node[ElementNodeId_s[ie][3]][j] + Disp_grad_1_on_Node[ElementNodeId_s[ie][7]][j]);
			Disp_grad_1_on_Node[ElementNodeId_s[ie][16]][j] = 0.5*(Disp_grad_1_on_Node[ElementNodeId_s[ie][4]][j] + Disp_grad_1_on_Node[ElementNodeId_s[ie][5]][j]);
			Disp_grad_1_on_Node[ElementNodeId_s[ie][17]][j] = 0.5*(Disp_grad_1_on_Node[ElementNodeId_s[ie][5]][j] + Disp_grad_1_on_Node[ElementNodeId_s[ie][6]][j]);
			Disp_grad_1_on_Node[ElementNodeId_s[ie][18]][j] = 0.5*(Disp_grad_1_on_Node[ElementNodeId_s[ie][6]][j] + Disp_grad_1_on_Node[ElementNodeId_s[ie][7]][j]);
			Disp_grad_1_on_Node[ElementNodeId_s[ie][19]][j] = 0.5*(Disp_grad_1_on_Node[ElementNodeId_s[ie][7]][j] + Disp_grad_1_on_Node[ElementNodeId_s[ie][4]][j]);

			Disp_grad_2_on_Node[ElementNodeId_s[ie][8]][j]  = 0.5*(Disp_grad_2_on_Node[ElementNodeId_s[ie][0]][j] + Disp_grad_2_on_Node[ElementNodeId_s[ie][1]][j]);
			Disp_grad_2_on_Node[ElementNodeId_s[ie][9]][j]  = 0.5*(Disp_grad_2_on_Node[ElementNodeId_s[ie][1]][j] + Disp_grad_2_on_Node[ElementNodeId_s[ie][2]][j]);
			Disp_grad_2_on_Node[ElementNodeId_s[ie][10]][j] = 0.5*(Disp_grad_2_on_Node[ElementNodeId_s[ie][2]][j] + Disp_grad_2_on_Node[ElementNodeId_s[ie][3]][j]);
			Disp_grad_2_on_Node[ElementNodeId_s[ie][11]][j] = 0.5*(Disp_grad_2_on_Node[ElementNodeId_s[ie][3]][j] + Disp_grad_2_on_Node[ElementNodeId_s[ie][0]][j]);
			Disp_grad_2_on_Node[ElementNodeId_s[ie][12]][j] = 0.5*(Disp_grad_2_on_Node[ElementNodeId_s[ie][0]][j] + Disp_grad_2_on_Node[ElementNodeId_s[ie][4]][j]);
			Disp_grad_2_on_Node[ElementNodeId_s[ie][13]][j] = 0.5*(Disp_grad_2_on_Node[ElementNodeId_s[ie][1]][j] + Disp_grad_2_on_Node[ElementNodeId_s[ie][5]][j]);
			Disp_grad_2_on_Node[ElementNodeId_s[ie][14]][j] = 0.5*(Disp_grad_2_on_Node[ElementNodeId_s[ie][2]][j] + Disp_grad_2_on_Node[ElementNodeId_s[ie][6]][j]);
			Disp_grad_2_on_Node[ElementNodeId_s[ie][15]][j] = 0.5*(Disp_grad_2_on_Node[ElementNodeId_s[ie][3]][j] + Disp_grad_2_on_Node[ElementNodeId_s[ie][7]][j]);
			Disp_grad_2_on_Node[ElementNodeId_s[ie][16]][j] = 0.5*(Disp_grad_2_on_Node[ElementNodeId_s[ie][4]][j] + Disp_grad_2_on_Node[ElementNodeId_s[ie][5]][j]);
			Disp_grad_2_on_Node[ElementNodeId_s[ie][17]][j] = 0.5*(Disp_grad_2_on_Node[ElementNodeId_s[ie][5]][j] + Disp_grad_2_on_Node[ElementNodeId_s[ie][6]][j]);
			Disp_grad_2_on_Node[ElementNodeId_s[ie][18]][j] = 0.5*(Disp_grad_2_on_Node[ElementNodeId_s[ie][6]][j] + Disp_grad_2_on_Node[ElementNodeId_s[ie][7]][j]);
			Disp_grad_2_on_Node[ElementNodeId_s[ie][19]][j] = 0.5*(Disp_grad_2_on_Node[ElementNodeId_s[ie][7]][j] + Disp_grad_2_on_Node[ElementNodeId_s[ie][4]][j]);


			Disp_grad_3_on_Node[ElementNodeId_s[ie][8]][j]  = 0.5*(Disp_grad_3_on_Node[ElementNodeId_s[ie][0]][j] + Disp_grad_3_on_Node[ElementNodeId_s[ie][1]][j]);
			Disp_grad_3_on_Node[ElementNodeId_s[ie][9]][j]  = 0.5*(Disp_grad_3_on_Node[ElementNodeId_s[ie][1]][j] + Disp_grad_3_on_Node[ElementNodeId_s[ie][2]][j]);
			Disp_grad_3_on_Node[ElementNodeId_s[ie][10]][j] = 0.5*(Disp_grad_3_on_Node[ElementNodeId_s[ie][2]][j] + Disp_grad_3_on_Node[ElementNodeId_s[ie][3]][j]);
			Disp_grad_3_on_Node[ElementNodeId_s[ie][11]][j] = 0.5*(Disp_grad_3_on_Node[ElementNodeId_s[ie][3]][j] + Disp_grad_3_on_Node[ElementNodeId_s[ie][0]][j]);
			Disp_grad_3_on_Node[ElementNodeId_s[ie][12]][j] = 0.5*(Disp_grad_3_on_Node[ElementNodeId_s[ie][0]][j] + Disp_grad_3_on_Node[ElementNodeId_s[ie][4]][j]);
			Disp_grad_3_on_Node[ElementNodeId_s[ie][13]][j] = 0.5*(Disp_grad_3_on_Node[ElementNodeId_s[ie][1]][j] + Disp_grad_3_on_Node[ElementNodeId_s[ie][5]][j]);
			Disp_grad_3_on_Node[ElementNodeId_s[ie][14]][j] = 0.5*(Disp_grad_3_on_Node[ElementNodeId_s[ie][2]][j] + Disp_grad_3_on_Node[ElementNodeId_s[ie][6]][j]);
			Disp_grad_3_on_Node[ElementNodeId_s[ie][15]][j] = 0.5*(Disp_grad_3_on_Node[ElementNodeId_s[ie][3]][j] + Disp_grad_3_on_Node[ElementNodeId_s[ie][7]][j]);
			Disp_grad_3_on_Node[ElementNodeId_s[ie][16]][j] = 0.5*(Disp_grad_3_on_Node[ElementNodeId_s[ie][4]][j] + Disp_grad_3_on_Node[ElementNodeId_s[ie][5]][j]);
			Disp_grad_3_on_Node[ElementNodeId_s[ie][17]][j] = 0.5*(Disp_grad_3_on_Node[ElementNodeId_s[ie][5]][j] + Disp_grad_3_on_Node[ElementNodeId_s[ie][6]][j]);
			Disp_grad_3_on_Node[ElementNodeId_s[ie][18]][j] = 0.5*(Disp_grad_3_on_Node[ElementNodeId_s[ie][6]][j] + Disp_grad_3_on_Node[ElementNodeId_s[ie][7]][j]);
			Disp_grad_3_on_Node[ElementNodeId_s[ie][19]][j] = 0.5*(Disp_grad_3_on_Node[ElementNodeId_s[ie][7]][j] + Disp_grad_3_on_Node[ElementNodeId_s[ie][4]][j]);

			W_1_on_Node[ElementNodeId_s[ie][8]][j]  = 0.5*(W_1_on_Node[ElementNodeId_s[ie][0]][j] + W_1_on_Node[ElementNodeId_s[ie][1]][j]);
			W_1_on_Node[ElementNodeId_s[ie][9]][j]  = 0.5*(W_1_on_Node[ElementNodeId_s[ie][1]][j] + W_1_on_Node[ElementNodeId_s[ie][2]][j]);
			W_1_on_Node[ElementNodeId_s[ie][10]][j] = 0.5*(W_1_on_Node[ElementNodeId_s[ie][2]][j] + W_1_on_Node[ElementNodeId_s[ie][3]][j]);
			W_1_on_Node[ElementNodeId_s[ie][11]][j] = 0.5*(W_1_on_Node[ElementNodeId_s[ie][3]][j] + W_1_on_Node[ElementNodeId_s[ie][0]][j]);
			W_1_on_Node[ElementNodeId_s[ie][12]][j] = 0.5*(W_1_on_Node[ElementNodeId_s[ie][0]][j] + W_1_on_Node[ElementNodeId_s[ie][4]][j]);
			W_1_on_Node[ElementNodeId_s[ie][13]][j] = 0.5*(W_1_on_Node[ElementNodeId_s[ie][1]][j] + W_1_on_Node[ElementNodeId_s[ie][5]][j]);
			W_1_on_Node[ElementNodeId_s[ie][14]][j] = 0.5*(W_1_on_Node[ElementNodeId_s[ie][2]][j] + W_1_on_Node[ElementNodeId_s[ie][6]][j]);
			W_1_on_Node[ElementNodeId_s[ie][15]][j] = 0.5*(W_1_on_Node[ElementNodeId_s[ie][3]][j] + W_1_on_Node[ElementNodeId_s[ie][7]][j]);
			W_1_on_Node[ElementNodeId_s[ie][16]][j] = 0.5*(W_1_on_Node[ElementNodeId_s[ie][4]][j] + W_1_on_Node[ElementNodeId_s[ie][5]][j]);
			W_1_on_Node[ElementNodeId_s[ie][17]][j] = 0.5*(W_1_on_Node[ElementNodeId_s[ie][5]][j] + W_1_on_Node[ElementNodeId_s[ie][6]][j]);
			W_1_on_Node[ElementNodeId_s[ie][18]][j] = 0.5*(W_1_on_Node[ElementNodeId_s[ie][6]][j] + W_1_on_Node[ElementNodeId_s[ie][7]][j]);
			W_1_on_Node[ElementNodeId_s[ie][19]][j] = 0.5*(W_1_on_Node[ElementNodeId_s[ie][7]][j] + W_1_on_Node[ElementNodeId_s[ie][4]][j]);

			W_2_on_Node[ElementNodeId_s[ie][8]][j]  = 0.5*(W_2_on_Node[ElementNodeId_s[ie][0]][j] + W_2_on_Node[ElementNodeId_s[ie][1]][j]);
			W_2_on_Node[ElementNodeId_s[ie][9]][j]  = 0.5*(W_2_on_Node[ElementNodeId_s[ie][1]][j] + W_2_on_Node[ElementNodeId_s[ie][2]][j]);
			W_2_on_Node[ElementNodeId_s[ie][10]][j] = 0.5*(W_2_on_Node[ElementNodeId_s[ie][2]][j] + W_2_on_Node[ElementNodeId_s[ie][3]][j]);
			W_2_on_Node[ElementNodeId_s[ie][11]][j] = 0.5*(W_2_on_Node[ElementNodeId_s[ie][3]][j] + W_2_on_Node[ElementNodeId_s[ie][0]][j]);
			W_2_on_Node[ElementNodeId_s[ie][12]][j] = 0.5*(W_2_on_Node[ElementNodeId_s[ie][0]][j] + W_2_on_Node[ElementNodeId_s[ie][4]][j]);
			W_2_on_Node[ElementNodeId_s[ie][13]][j] = 0.5*(W_2_on_Node[ElementNodeId_s[ie][1]][j] + W_2_on_Node[ElementNodeId_s[ie][5]][j]);
			W_2_on_Node[ElementNodeId_s[ie][14]][j] = 0.5*(W_2_on_Node[ElementNodeId_s[ie][2]][j] + W_2_on_Node[ElementNodeId_s[ie][6]][j]);
			W_2_on_Node[ElementNodeId_s[ie][15]][j] = 0.5*(W_2_on_Node[ElementNodeId_s[ie][3]][j] + W_2_on_Node[ElementNodeId_s[ie][7]][j]);
			W_2_on_Node[ElementNodeId_s[ie][16]][j] = 0.5*(W_2_on_Node[ElementNodeId_s[ie][4]][j] + W_2_on_Node[ElementNodeId_s[ie][5]][j]);
			W_2_on_Node[ElementNodeId_s[ie][17]][j] = 0.5*(W_2_on_Node[ElementNodeId_s[ie][5]][j] + W_2_on_Node[ElementNodeId_s[ie][6]][j]);
			W_2_on_Node[ElementNodeId_s[ie][18]][j] = 0.5*(W_2_on_Node[ElementNodeId_s[ie][6]][j] + W_2_on_Node[ElementNodeId_s[ie][7]][j]);
			W_2_on_Node[ElementNodeId_s[ie][19]][j] = 0.5*(W_2_on_Node[ElementNodeId_s[ie][7]][j] + W_2_on_Node[ElementNodeId_s[ie][4]][j]);

			W_3_on_Node[ElementNodeId_s[ie][8]][j]  = 0.5*(W_3_on_Node[ElementNodeId_s[ie][0]][j] + W_3_on_Node[ElementNodeId_s[ie][1]][j]);
			W_3_on_Node[ElementNodeId_s[ie][9]][j]  = 0.5*(W_3_on_Node[ElementNodeId_s[ie][1]][j] + W_3_on_Node[ElementNodeId_s[ie][2]][j]);
			W_3_on_Node[ElementNodeId_s[ie][10]][j] = 0.5*(W_3_on_Node[ElementNodeId_s[ie][2]][j] + W_3_on_Node[ElementNodeId_s[ie][3]][j]);
			W_3_on_Node[ElementNodeId_s[ie][11]][j] = 0.5*(W_3_on_Node[ElementNodeId_s[ie][3]][j] + W_3_on_Node[ElementNodeId_s[ie][0]][j]);
			W_3_on_Node[ElementNodeId_s[ie][12]][j] = 0.5*(W_3_on_Node[ElementNodeId_s[ie][0]][j] + W_3_on_Node[ElementNodeId_s[ie][4]][j]);
			W_3_on_Node[ElementNodeId_s[ie][13]][j] = 0.5*(W_3_on_Node[ElementNodeId_s[ie][1]][j] + W_3_on_Node[ElementNodeId_s[ie][5]][j]);
			W_3_on_Node[ElementNodeId_s[ie][14]][j] = 0.5*(W_3_on_Node[ElementNodeId_s[ie][2]][j] + W_3_on_Node[ElementNodeId_s[ie][6]][j]);
			W_3_on_Node[ElementNodeId_s[ie][15]][j] = 0.5*(W_3_on_Node[ElementNodeId_s[ie][3]][j] + W_3_on_Node[ElementNodeId_s[ie][7]][j]);
			W_3_on_Node[ElementNodeId_s[ie][16]][j] = 0.5*(W_3_on_Node[ElementNodeId_s[ie][4]][j] + W_3_on_Node[ElementNodeId_s[ie][5]][j]);
			W_3_on_Node[ElementNodeId_s[ie][17]][j] = 0.5*(W_3_on_Node[ElementNodeId_s[ie][5]][j] + W_3_on_Node[ElementNodeId_s[ie][6]][j]);
			W_3_on_Node[ElementNodeId_s[ie][18]][j] = 0.5*(W_3_on_Node[ElementNodeId_s[ie][6]][j] + W_3_on_Node[ElementNodeId_s[ie][7]][j]);
			W_3_on_Node[ElementNodeId_s[ie][19]][j] = 0.5*(W_3_on_Node[ElementNodeId_s[ie][7]][j] + W_3_on_Node[ElementNodeId_s[ie][4]][j]);	   	}
  	}

  	for(ip=0; ip<NNodes; ip++){
  		if( isnan(JW_on_Node[ip][0]) != 0 ){
	  		//printf("%d %d\n",ip, N_SP);
	  		//printf("PCID %d %e %e %e\n", ip, Coord_of_PatchCenter[ip][0], Coord_of_PatchCenter[ip][1], Coord_of_PatchCenter[ip][2]);
	  		//printf(" %e %e %e", Coord_of_PatchCenter[ip][0], Coord_of_PatchCenter[ip][1], Coord_of_PatchCenter[ip][2]);
	  		for(i=0; i<1; i++){
				//printf(" %e", JW_on_Node[ip][i]);
	  		}
	  		//printf("\n");
		}
  	}
}

/*****************************************************************************
   			積分範囲のリスト作成
*******************************************************************************/
void Make_J_Integral_List_T10(){
	int ie,ip,icp,i,j,D_check;
	double xc,yc,zc;
  double xxc,yyc,zzc;
  double d_x,d_y,d_z,D;

 ////////要素リスト作成//////////
	for(ie=0;ie<NElements;ie++){
		D_check = 0;
		xc=0;
	  yc=0;
		zc=0;
		//要素の重心点
		for(j=0;j<4;j++){
			xc+=(Node_Coordinate[ElementNodeId_s[ie][j]][0])/4;
			yc+=(Node_Coordinate[ElementNodeId_s[ie][j]][1])/4;
			zc+=(Node_Coordinate[ElementNodeId_s[ie][j]][2])/4;
		}

		for(icp=0;icp<CrackNodes_DI;icp++){
			xxc = Node_Coordinate[List_CrackNode_DI[icp]][0];
			yyc = Node_Coordinate[List_CrackNode_DI[icp]][1];
			zzc = Node_Coordinate[List_CrackNode_DI[icp]][2];
			d_x=xxc-xc;
    	d_y=yyc-yc;
    	d_z=zzc-zc;
			D=sqrt(d_x*d_x + d_y*d_y + d_z*d_z);
			if(D <= RRR) {
				D_check ++;
				break;
			}
		}
		if(D_check >= 1){
			List_Element_DI[NElements_DI] = ie;
			NElements_DI ++;
		}
	}

 ////////節点リスト作成//////////
	for(ip=0;ip<NNodes;ip++){
		D_check = 0;
		xc=Node_Coordinate[ip][0];
		yc=Node_Coordinate[ip][1];
		zc=Node_Coordinate[ip][2];

		for(icp=0;icp<CrackNodes_DI;icp++){
			xxc = Node_Coordinate[List_CrackNode_DI[icp]][0];
			yyc = Node_Coordinate[List_CrackNode_DI[icp]][1];
			zzc = Node_Coordinate[List_CrackNode_DI[icp]][2];
			d_x=xxc-xc;
    	d_y=yyc-yc;
    	d_z=zzc-zc;
			D=sqrt(d_x*d_x + d_y*d_y + d_z*d_z);
			if(D <= RRR) {
				D_check ++;
				break;
			}
		}
		if(D_check >= 1){
			List_Node_DI[NNodes_DI]=ip;
			NNodes_DI ++;
		}
	}

	printf("\n----------------------\nList_Element_DI = %d\nList_Node_DI = %d\n----------------------\n\n",NElements_DI,NNodes_DI);

}

//ElementList可視化用
void Oudput_Integration_Elements_List_T10(){
	FILE *fp;
  	fp = fopen("Integration_Elements_List_Pre.out", "w");
  	int ie;

  	fprintf(fp,"%d\n",NElements_DI);
  	for(ie=0;ie<NElements_DI;ie++){
  		fprintf(fp,"%d\n",List_Element_DI[ie]);
  	}
  	fclose(fp);
}

//int Make_Value_on_Node(int data)
int Make_Value_on_Node(int data){

  	int N_COMPONENT=0;

	switch(data){
		case 1:
			N_COMPONENT=1;
			printf("Make_Value_on_Node(W)\n");
			break;
		case 2:
			N_COMPONENT=9;
			printf("Make_Value_on_Node(Pai)\n");
			break;
		case 3:
			N_COMPONENT=3;
			printf("Make_Value_on_Node(Disp_grad_1)\n");
			break;
		case 4:
			N_COMPONENT=3;
			printf("Make_Value_on_Node(Disp_grad_2)\n");
			break;
		case 5:
			N_COMPONENT=3;
			printf("Make_Value_on_Node(Disp_grad_3)\n");
			break;
		case 6:
			N_COMPONENT=3;
			printf("Make_Value_on_Node(W_1)\n");
			break;
		case 7:
			N_COMPONENT=3;
			printf("Make_Value_on_Node(W_2)\n");
			break;
		case 8:
			N_COMPONENT=3;
			printf("Make_Value_on_Node(W_3)\n");
			break;
	}

  	int ie, ig, ip, ipl, id, in, ic, i,j,k,n;
  	int N_SP;

  	if(method==1){//local least square method linear
  		for(ipl=0; ipl<NNodes_DI; ipl++){
  			ip = List_Node_DI[ipl];
			  SP_all_Gausspoint_in_eachelement_linear(data, N_COMPONENT, ip);//節点(PC)を含む各要素の全ガウス点を参照点として計算
	  		if(ipl!=0)
	  			if ((ipl/DISPLAY)*DISPLAY == ipl){
					printf ("%d Nodes End\n", ipl);
	  			}
  		}
  		if(Number_of_Nodes_on_Elements==10){
  			printf("make intermediate node data\n");
   			data_to_intermediate_node_T10();//頂点節点に与えたデータから中間節点にデータを与える。中間節点は頂点節点の平均。
  		}
  		else if(Number_of_Nodes_on_Elements==20){
  			printf("make intermediate node data\n");
   			data_to_intermediate_node_H20();//頂点節点に与えたデータから中間節点にデータを与える。中間節点は頂点節点の平均。
  		}
  		printf("local least square method linear finished\n");
   			
  	}

  	else if(method==2){//local least square method quad
  		for(ipl=0; ipl<NNodes_DI; ipl++){
  			ip = List_Node_DI[ipl];
			  SP_all_Gausspoint_in_eachelement_quad(data, N_COMPONENT, ip);//節点(PC)を含む各要素の全ガウス点を参照点として計算
  			//printf("ip=%d\n", ip);
	  		if(ipl!=0){
	  			if ((ipl/DISPLAY)*DISPLAY == ipl){
					printf ("%d Nodes End\n", ipl);
	  			}
	  		}
  		}
  		printf("local least square method quad finished\n");
  	}

  	else if(method==3){//average
  		for(ipl=0; ipl<NNodes_DI; ipl++){
  			ip = List_Node_DI[ipl];
   			SP_Average_one_Gausspoint_in_eachelement(data, N_COMPONENT, ip);//節点(PC)を含む各要素の節点に一番近いガウス点の平均値を節点に与える
  			//printf("ip=%d\n", ip);
	  		if(ipl!=0)
	  			if ((ipl/DISPLAY)*DISPLAY == ipl){
					printf ("%d Nodes End\n", ipl);
	  			}
  		}
  		if(Number_of_Nodes_on_Elements==10){
  			printf("make intermediate node data\n");
   			data_to_intermediate_node_T10();//頂点節点に与えたデータから中間節点にデータを与える。中間節点は頂点節点の平均。
  		}
  		else if(Number_of_Nodes_on_Elements==20){
  			printf("make intermediate node data\n");
   			data_to_intermediate_node_H20();//頂点節点に与えたデータから中間節点にデータを与える。中間節点は頂点節点の平均。
  		}
  	}

  	output_node_data(data, N_COMPONENT);
  	//output_node_data(data, N_COMPONENT, fileName);
}

/**************************************************************
                     	 データ出力
*************************************************************/
/***********************************************************************************
データ出力フォーマットはすべて同じ書き方。
(ex)
221 (←総接点数)
2.98988033533469703684e+01 4.37624549383086502985e+00 2.40772688113523081199e+00 
4.37624549383086502985e+00 1.37860520342081940726e+02 8.12191937326521440355e+00 
2.40772688113523081199e+00 8.12191937326521440355e+00 4.35166979931521922254e+01 
(↑節点番号0の1stPKの全9成分)
2.98988033533450483503e+01 4.37624549383060923446e+00 -2.40772688113500121787e+00 
4.37624549383060923446e+00 1.37860520342073414213e+02 -8.12191937326664969987e+00 
-2.40772688113500121787e+00 -8.12191937326664969987e+00 4.35166979931494779521e+01 
(↑節点番号1の1stPKの全9成分)
***********************************************************************************/
int output_node_data(int data, int N_COMPONENT){

  	int i,j;
  	int temp_int;
  	FILE *fp;

	if(data==1){
		fp=fopen("W_on_Node.dat", "w");
		fprintf(fp, "%d\n", NNodes);

		for(i=0; i<NNodes; i++){
			for(j=0; j<N_COMPONENT; j++){
				fprintf(fp,"%21.20le ", JW_on_Node[i][j]);
			}
			fprintf(fp,"\n");
		}
		if(method==1)		fprintf(fp,"local least square method linear\n");
		else if(method==2)	fprintf(fp,"local least square method quad\n");
		else if(method==3)	fprintf(fp,"average\n");
		if(deform==1)		fprintf(fp,"small deformation\n");
		else if(deform==2)	fprintf(fp,"large deformation\n");
		fclose(fp);
	}

	if(data==2){
		fp=fopen("Pai_on_Node.dat", "w");
		fprintf(fp, "%d\n", NNodes);
		for(i=0; i<NNodes; i++){
			for(j=0; j<N_COMPONENT; j++){
				fprintf(fp,"%21.20le ", Pai_on_Node[i][j]);
				if(j%3==2)fprintf(fp,"\n");
			}
			fprintf(fp,"\n");
		}
		if(method==1)		fprintf(fp,"local least square method linear\n");
		else if(method==2)	fprintf(fp,"local least square method quad\n");
		else if(method==3)	fprintf(fp,"average\n");
		if(deform==1)		fprintf(fp,"small deformation\n");
		else if(deform==2)	fprintf(fp,"large deformation\n");
		fclose(fp);
	}

	if(data==3){
		fp=fopen("Disp_grad_1_on_Node.dat", "w");
		fprintf(fp, "%d\n", NNodes);
		for(i=0; i<NNodes; i++){
			for(j=0; j<N_COMPONENT; j++){
				fprintf(fp,"%21.20le ", Disp_grad_1_on_Node[i][j]);
			}
			fprintf(fp,"\n");
		}
		if(method==1)		fprintf(fp,"local least square method linear\n");
		else if(method==2)	fprintf(fp,"local least square method quad\n");
		else if(method==3)	fprintf(fp,"average\n");
		if(deform==1)		fprintf(fp,"small deformation\n");
		else if(deform==2)	fprintf(fp,"large deformation\n");
		fclose(fp);
	}

	if(data==4){
		fp=fopen("Disp_grad_2_on_Node.dat", "w");
		fprintf(fp, "%d\n", NNodes);
		for(i=0; i<NNodes; i++){
			for(j=0; j<N_COMPONENT; j++){
				fprintf(fp,"%21.20le ", Disp_grad_2_on_Node[i][j]);
			}
			fprintf(fp,"\n");
		}
		if(method==1)		fprintf(fp,"local least square method linear\n");
		else if(method==2)	fprintf(fp,"local least square method quad\n");
		else if(method==3)	fprintf(fp,"average\n");
		if(deform==1)		fprintf(fp,"small deformation\n");
		else if(deform==2)	fprintf(fp,"large deformation\n");
		fclose(fp);
	}

	if(data==5){
		fp=fopen("Disp_grad_3_on_Node.dat", "w");
		fprintf(fp, "%d\n", NNodes);
		for(i=0; i<NNodes; i++){
			for(j=0; j<N_COMPONENT; j++){
				fprintf(fp,"%21.20le ", Disp_grad_3_on_Node[i][j]);
			}
			fprintf(fp,"\n");
		}
		if(method==1)		fprintf(fp,"local least square method linear\n");
		else if(method==2)	fprintf(fp,"local least square method quad\n");
		else if(method==3)	fprintf(fp,"average\n");
		if(deform==1)		fprintf(fp,"small deformation\n");
		else if(deform==2)	fprintf(fp,"large deformation\n");
		fclose(fp);
	}

	if(data==6){
		fp=fopen("W_1_on_Node.dat", "w");
		fprintf(fp, "%d\n", NNodes);
		for(i=0; i<NNodes; i++){
			for(j=0; j<N_COMPONENT; j++){
				fprintf(fp,"%21.20le ", W_1_on_Node[i][j]);
			}
			fprintf(fp,"\n");
		}
		if(method==1)		fprintf(fp,"local least square method linear\n");
		else if(method==2)	fprintf(fp,"local least square method quad\n");
		else if(method==3)	fprintf(fp,"average\n");
		if(deform==1)		fprintf(fp,"small deformation\n");
		else if(deform==2)	fprintf(fp,"large deformation\n");
		fclose(fp);
	}

	if(data==7){
		fp=fopen("W_2_on_Node.dat", "w");
		fprintf(fp, "%d\n", NNodes);
		for(i=0; i<NNodes; i++){
			for(j=0; j<N_COMPONENT; j++){
				fprintf(fp,"%21.20le ", W_2_on_Node[i][j]);
			}
			fprintf(fp,"\n");
		}
		if(method==1)		fprintf(fp,"local least square method linear\n");
		else if(method==2)	fprintf(fp,"local least square method quad\n");
		else if(method==3)	fprintf(fp,"average\n");
		if(deform==1)		fprintf(fp,"small deformation\n");
		else if(deform==2)	fprintf(fp,"large deformation\n");
		fclose(fp);
	}

	if(data==8){
		fp=fopen("W_3_on_Node.dat", "w");
		fprintf(fp, "%d\n", NNodes);
		for(i=0; i<NNodes; i++){
			for(j=0; j<N_COMPONENT; j++){
				fprintf(fp,"%21.20le ", W_3_on_Node[i][j]);
				}
			fprintf(fp,"\n");
		}
		if(method==1)		fprintf(fp,"local least square method linear\n");
		else if(method==2)	fprintf(fp,"local least square method quad\n");
		else if(method==3)	fprintf(fp,"average\n");
		if(deform==1)		fprintf(fp,"small deformation\n");
		else if(deform==2)	fprintf(fp,"large deformation\n");
		fclose(fp);
	}
}

/**************************************************************
                     	 J Integral main
*************************************************************/
void Make_data_J_Integral(){

	switch(Number_of_Nodes_on_Elements){
	  case 4:
	  	printf("Tetra 4 mode\n");
	  	Make_Disp_grad_2_Tetra_4();
	  	Make_1st_PK_stress();
	  	make_Gxi_coord_Tetra_4();
	  	Make_JW();
	  	Make_EMT();
	  	break;
	  case 10:
	  	printf("Tetra 10 mode\n");
			Make_J_Integral_List_T10();
			Oudput_Integration_Elements_List_T10();
	  	Make_Disp_grad_2_Tetra_10();
	  	Make_1st_PK_stress();
	  	make_Gxi_coord_Tetra_10();
	  	Make_JW();
	  	Make_EMT();
	  	break;
	  case 8:
	  	printf("Hexa 8 mode\n");
	  	Make_Disp_grad_2_Hexa_8();
	  	Make_1st_PK_stress();
	  	make_Gxi_coord_Hexa_8();
		Make_JW();
		Make_EMT();
	  	break;
	  case 20:
	  	printf("Hexa 20 mode\n");
	  	Make_Disp_grad_2_Hexa_20();
	  	Make_1st_PK_stress();
	  	make_Gxi_coord_Hexa_20();
		Make_JW();
		Make_EMT();
	  	break;
	}		

	if(method==1){
		Make_Value_on_Node(1);//(W)
		Make_Value_on_Node(2);//(Pai)
		Make_Value_on_Node(3);//(Disp_grad_1)
		Make_Value_on_Node(4);//(Disp_grad_2)
		Make_Value_on_Node(5);//(Disp_grad_3)
		//Make_Value_on_Node(6);//(W_1) //ADVENTURE_J-Integral_largedef_ver5_2.xでsecond_eq_type=3を用いる時はコメントアウトを外す
		//Make_Value_on_Node(7);//(W_2) //ADVENTURE_J-Integral_largedef_ver5_2.xでsecond_eq_type=3を用いる時はコメントアウトを外す
		//Make_Value_on_Node(8);//(W_3) //ADVENTURE_J-Integral_largedef_ver5_2.xでsecond_eq_type=3を用いる時はコメントアウトを外す
	}

	if(method==2){
		Make_Value_on_Node(1);//(W)
		////Make_Value_on_Node(2);//(Pai)  局所最小二乗法の二次基底では使用しない
		////Make_Value_on_Node(3);//(Disp_grad_1)局所最小二乗法の二次基底では使用しない
		////Make_Value_on_Node(4);//(Disp_grad_2)局所最小二乗法の二次基底では使用しない
		////Make_Value_on_Node(5);//(Disp_grad_3)局所最小二乗法の二次基底では使用しない
		//Make_Value_on_Node(6);//(W_1) //ADVENTURE_J-Integral_largedef_ver5_2.xでsecond_eq_type=3を用いる時はコメントアウトを外す
		//Make_Value_on_Node(7);//(W_2) //ADVENTURE_J-Integral_largedef_ver5_2.xでsecond_eq_type=3を用いる時はコメントアウトを外す
		//Make_Value_on_Node(8);//(W_3) //ADVENTURE_J-Integral_largedef_ver5_2.xでsecond_eq_type=3を用いる時はコメントアウトを外す
	}

	if(method==3){
		Make_Value_on_Node(1);//(W)
		Make_Value_on_Node(2);//(Pai)
		Make_Value_on_Node(3);//(Disp_grad_1)
		Make_Value_on_Node(4);//(Disp_grad_2)
		Make_Value_on_Node(5);//(Disp_grad_3)
		//Make_Value_on_Node(6);//(W_1) //ADVENTURE_J-Integral_largedef_ver5_2.xでsecond_eq_type=3を用いる時はコメントアウトを外す
		//Make_Value_on_Node(7);//(W_2) //ADVENTURE_J-Integral_largedef_ver5_2.xでsecond_eq_type=3を用いる時はコメントアウトを外す
		//Make_Value_on_Node(8);//(W_3) //ADVENTURE_J-Integral_largedef_ver5_2.xでsecond_eq_type=3を用いる時はコメントアウトを外す
	}

}

/**************************************************************
                     	 データ出力系     旧J積分コードとガウス点の順番を合わせて出力してる（convertNとか）
*************************************************************/
int Output_Stress_Data(){

  	int i,j,k,l;
	FILE *fp;
	fp = fopen("test_Stress", "w");
	int convertN[4]={2,3,0,1},kk;

	fprintf(fp,"%d\n", NElements);
	//応力
	for( i = 0; i < NElements; i++ ){
		for( kk = 0; kk < N_IntegrationPoint; kk++){
			if(Number_of_Nodes_on_Elements==10) k=convertN[kk];
			else k=kk;
		  	for( j = 0; j < DIMENSION; j++ ){
				for( l = 0; l < DIMENSION; l++ ){
			  		fprintf(fp,"%5.20le ",Stress[i][k][j][l]);
			  		if(l==2) fprintf(fp,"\n");
				}
			}
			fprintf(fp,"\n");
		}
	}
  	fclose(fp);
}

int Output_Disp_grad_Data(){

  	int i,j,k,l;
  	int N,e;
	FILE *fp;
	fp = fopen("Disp_grad", "w");
	int convertN[4]={2,3,0,1},NN;

	fprintf(fp,"%d\n", NElements); 
	for( e = 0; e < NElements; e++ ){
		for( NN = 0; NN < N_IntegrationPoint; NN++ ){
			if(Number_of_Nodes_on_Elements==10) N=convertN[NN];
			else N=NN;

		  	for( i = 0; i < DIMENSION; i++ ){
			  	fprintf(fp,"element %d iGp:%d | %5.20le %5.20le %5.20le\n",e,N, Disp_grad_1[e][N][i], Disp_grad_2[e][N][i], Disp_grad_3[e][N][i]);
 				if(i==2) fprintf(fp,"\n");
		  	}
		}
	}
  	fclose(fp);
}

int Output_Def_grad_Data(){

  	int i,j,k,l;
 	int N,e;
  	FILE *fp;
  	fp = fopen("test_Def_grad", "w");
  	int convertN[4]={2,3,0,1},NN;

	fprintf(fp,"%d\n", NElements);
	for( e = 0; e < NElements; e++ ){
		for( NN = 0; NN < N_IntegrationPoint; NN++ ){
			if(Number_of_Nodes_on_Elements==10)N=convertN[NN];
			else N=NN;
	  		for( i = 0; i < DIMENSION; i++ ){
				for( j = 0; j < DIMENSION; j++ ){
	  				fprintf(fp,"%5.20le ", Def_grad[e][N][i][j]);
	  				if(j==2) fprintf(fp,"\n");
				}
				if(i==2) fprintf(fp,"\n");
	  		}
		}
  	}
  	fclose(fp);
}

int Output_Def_grad_In_Data(){

  	int i,j,k,l;
  	int N,e;
  	FILE *fp;
  	fp = fopen("test_Def_grad_In", "w");
  	int convertN[4]={2,3,0,1},NN;

	fprintf(fp,"%d\n", NElements);
  	for( e = 0; e < NElements; e++ ){
		for( NN = 0; NN < N_IntegrationPoint; NN++ ){
			if(Number_of_Nodes_on_Elements==10)N=convertN[NN];
			else N=NN;
	  		for( i = 0; i < DIMENSION; i++ ){
				for( j = 0; j < DIMENSION; j++ ){
	  				fprintf(fp,"%5.20le ", Def_grad_In[e][N][i][j]);
	  				if(j==2) fprintf(fp,"\n");
				}
				if(i==2) fprintf(fp,"\n");
	  		}
		}
  	}
  	fclose(fp);
}

int Output_I(){

  	int i,j,k,l;
  	int N,e;
  	FILE *fp;
  	fp = fopen("test_I", "w");
  	int convertN[4]={2,3,0,1},NN;

	fprintf(fp,"%d", NElements);
  	for( e = 0; e < NElements; e++ ){
		for( NN = 0; NN < N_IntegrationPoint; NN++ ){
			if(Number_of_Nodes_on_Elements==10)N=convertN[NN];
			else N=NN;	
	  		for( i = 0; i < DIMENSION; i++ ){
				for( j = 0; j < DIMENSION; j++ ){
	  				fprintf(fp,"%5.20le ", I[e][N][i][j]);
	  				if(j==2) fprintf(fp,"\n");
				}
				if(i==2) fprintf(fp,"\n");
	  		}
		}
  	}
  fclose(fp);
}

int Output_Pai_Data(){

  	int i,j,k,l;
 	int N,e;
  	FILE *fp;
  	fp = fopen("test_Pai", "w");
  	int convertN[4]={2,3,0,1},NN;
	
	fprintf(fp,"%d\n", NElements);
  	for( e = 0; e < NElements; e++ ){
		for( NN = 0; NN < N_IntegrationPoint; NN++ ){
			if(Number_of_Nodes_on_Elements==10)N=convertN[NN];
			else N=NN;
	  		for( i = 0; i < DIMENSION; i++ ){
				for( j = 0; j < DIMENSION; j++ ){
	  				fprintf(fp,"%5.20le ", Pai[e][N][i][j]);
					//printf("PAI=%e\n",Pai[e][N][i][j]);
	  				if(j==2) fprintf(fp,"\n");
				}
				if(i==2) fprintf(fp,"\n");
	  		}
		}
  	}
  	fclose(fp);
}

/**************************************************************
                  			main
*************************************************************/
int main (int argc, char *argv[]) 
{

	if (argc != 10 ) {
			fprintf (stderr, 
			" Error\n"
			" argv[1]:.msh   argv[2]:.fgr\n"
			" argv[3]:ImputData from FEM\n"
			" argv[4]:0 or 1 (Select Displacement Mode)\n"
			" argv[5]:1 or 2 (Select Deformation Mode)\n"
			" argv[6]:1 or 2 or 3 (Select local square method or average)\n"
			" argv[7]:Crack Flag(upper or lower nodes)\n"
			" argv[8]:J_DI_input_data	(temp_J_int.dat)\n"
			" argv[9]:max RR \n");
			exit (1);
	}
	DispMode=atoi(argv[4]);//Xに変位量を足しこむかどうかの指定。足し込まない→0、足しこむ→1

	printf("-----------------------------------------------------------------\n[Selected Option]\n");
	deform = atoi(argv[5]);
	if(deform==1)printf("deform=%d small deformation\n",deform);
	else if(deform==2)printf("deform=%d large deformation\n",deform);
	else {printf("ERROR!! deform data argv[5] miss\n"); exit(1);}

	method = atoi(argv[6]);
	if(method==1)printf("method=%d local least square method linear\n",method);
	else if(method==2)printf("method=%d local least square method quad\n",method);
	else if(method==3)printf("method=%d average\n",method);
	else {printf("\nERROR!! method data argv[6] miss\n"); exit(1);}

	Max_RR = atof(argv[9]);
	printf("Max_RR=%lf\n",Max_RR);

	printf("-----------------------------------------------------------------\n[Mesh Property]\n");
	ReadFgr(argv[2]);	
	read_file_msh(argv[1]);
	reading_J_DI_input_data(argv[8]);

	printf("Number of Elements =%d\nNumber of Nodes=%d\n",NElements,NNodes);
	printf("-----------------------------------------------------------------\n");

	if(method==2){//use only local least square method quad
		read_Crack_flag(argv[7]);//局所最小二乗法を行うときにき裂をまたいで参照点をとらないようにするためのフラグ。節点ごとにある。
		make_Gauss_flag();//局所最小二乗法を行うときにき裂をまたいで参照点をとらないようにするためのフラグ。要素ごとにある。
	}

	printf("Reading FEM Result Data\n");
	Get_InputData_from_FEM(argv[3]);

	///////////////J_Integral//////////////////////////////////////////////
	printf("starting make data J integral\n");
	Make_data_J_Integral();
	printf("make data J-integral computation ended\n");
	Output_Stress_Data();
	Output_Pai_Data();
	Output_Def_grad_Data();
	Output_I();
}

