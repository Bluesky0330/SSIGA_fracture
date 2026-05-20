/*
2024/12/3
JintegralPointsはき裂前縁方向のパッチ数の等倍数であることを前提としている．
き裂前縁をη軸方向にそろえるモデリングを行う必要あり，縮退・特異パッチを使用する場合はξとζ方向を用いて極座標形式としている．
積分領域がパッチから出ないようdomainWidth < (JintegralPoints + 1) / 2  になるようにインプットを設定する必要あり．
特異パッチ，四辺形パッチともに要素を等分割していることを前提とした実装
特異パッチと四辺形パッチでq 関数の分布が異なる
J積分の実装は未完了
*/


// header
#include "_header.hpp"
#include "_IIM3D.hpp"


// input data の読みとり
void CrackDataSetting_3D::readCrackData()
{
    std::ifstream crackDataFile("J_int.txt");
    if (!crackDataFile) {
        std::cout << "Cannot open file: J_int.txt" << std::endl;
        exit(false);
    }

    // JintegralPoints:き裂前縁方向(η方向)のSIF計算点数
    crackDataFile >> domainWidth >> JintegralPoints;
    setJintegralPoints(JintegralPoints);

    std::cout << "Domain Width: " << domainWidth << std::endl;
    std::cout << "JintegralPoints: " << getJintegralPoints() << std::endl; 


    // き裂前縁(η)方向積分領域幅
    integralDomainNaturalWidth = 1.0 / domainWidth;                                     
    setIntegralDomainNaturalWidth(integralDomainNaturalWidth);
    cout << "IntegralDomainNaturalWidth: " << getIntegralDomainNaturalWidth() << endl;

    // き裂面の積分領域を設定 (0:上下両面, 1:どちらかのみ)
    crackDataFile >> model_size_para_3D;
    cout << "model_size_para_3D: " << getmodel_size_para_3D() << endl;

    // ξ方向の積分領域の大きさの設定
    crackDataFile >> domain_types;
    cout << "domain_types: " << getdomain_types() << endl;

    // ΔAの計算に置けるき裂前縁方向のき裂長さの計算に使用する領域分割数
    crackDataFile >> domain_area_divisions;
    cout << "domain_area_divisions: " << getdomain_area_divisions() << endl;

    // き裂先端パッチの種類(0: 縮退・特異パッチ, 1: 四辺形パッチ)
    int patchType;
    crackDataFile >> patchType;
    setpatchType(patchType);
    // cout << "patchType: " << patchType << endl;
    cout << "patchType: " << getpatchType() << endl;

    // read patch list
    int num_edge_crack_patch;  // き裂前縁方向のパッチ数
    int num_singular_patch;    // 特異パッチ数
    crackDataFile >> num_edge_crack_patch;
    crackDataFile >> num_singular_patch;
    setnum_edge_crack_patch(num_edge_crack_patch);
    setnum_singular_patch(num_singular_patch);
    // cout << "num_edge_crack_patch: " << num_edge_crack_patch << endl;
    cout << "num_edge_crack_patch: " << getnum_edge_crack_patch() << endl;
    // cout << "num_singular_patch: " << num_singular_patch << endl;
    cout << "num_singular_patch: " << getnum_singular_patch() << endl;

    // 特異パッチのIDを読み込む
    for (int i = 0; i < num_edge_crack_patch * num_singular_patch; i++) 
    {
        int patchID;
        crackDataFile >> patchID;
        patchList.push_back(patchID);
    }
    setpatchList(patchList);

    cout << "patch List: " << patchList.size() << endl;
    for (size_t i = 0; i < patchList.size(); i++) 
    {
        // cout << "patchList[" << i << "]: " << patchList[i] << endl;
        cout << "patch List[" << i << "]: " << getpatchList()[i] << endl;
    }
}

#if 1
// 計算パッチの各方向の要素数の取得
void CrackDataSetting_3D::setCrackTipPatchElement(information *info)
{   
    int local_init_patch = info->Total_Patch_to_mesh[Total_mesh - 1];  // ローカルパッチ番号の初期値
    int patchID = getpatchList()[0];
    patchID += local_init_patch;

    for (int dim = 0; dim < info->DIMENSION; dim  = dim + 2) 
    {
        vector<int> divisions(info->DIMENSION, 1);
        integralDomainDivisions[dim] = info->No_Control_point[patchID * info->DIMENSION + dim] - info->Order[patchID * info->DIMENSION + dim];
        // std::cout << "divisions[" << dim << "]: " << getIntegralDomainDivisions()[dim] << std::endl;
    }

    // for debug
    #if 1
    std::cout << "initial J-integral calculation patchID: " << patchID << std::endl;
    for (int dim = 0; dim < info->DIMENSION; dim++) 
    {
        std::cout << "divisions[" << dim << "]: " << getIntegralDomainDivisions()[dim] << std::endl;
    }
    #endif
}
#endif


void IIM_3D::computeIIM3D(CrackDataSetting_3D *crackData, information *info)
{
    // initialize
    vector<double> delta(info->DIMENSION * info->DIMENSION, 0.0);                    // クロネッカのデルタ
    for (int i = 0; i < info->DIMENSION; i++) 
        delta[i * info->DIMENSION + i] = 1.0;                                        // クロネッカのデルタの対角成分を1に設定


    // iterate over local patch elements
    for (int re = info->real_Total_Element_to_mesh[Total_mesh - 1]; re < info->real_Total_Element_to_mesh[Total_mesh]; re++)  // 現在参照パッチの要素数分の処理
    {
        int e = info->real_element[re];
        elementID_in_patch[info->Element_patch[e]].push_back(e);  // 各パッチに属する要素ID のリストを追加  
        ++elements_in_patch[info->Element_patch[e]];              // 各パッチに含まれる要素の数をカウント   
    }

    // for debug
    #if 0
    for (int re = info->real_Total_Element_to_mesh[Total_mesh - 1]; re < info->real_Total_Element_to_mesh[Total_mesh]; re++)  // 現在参照パッチの要素数分の処理
    {
        int e = info->real_element[re];
        printf("elementID_in_patch[%d]: %d\n", info->Element_patch[e], e);
        printf("elements_in_patch[%d]: %d\n", info->Element_patch[e], elements_in_patch[info->Element_patch[e]]);
    }
    #endif
    

    // compute J-integral at the crack front points
    std::unique_ptr<AuxiliaryField> auxField = std::make_unique<AuxiliaryField>();
    for (int Jcp = 0; Jcp < crackData->getJintegralPoints(); Jcp++)  
    {  
        // initialize J & M integral parameters
        for (int i = 0; i < crackData->getdomain_types(); i++) 
        {
            setJintegralScalar(Jcp * crackData->getdomain_types() + i, 0.0);
            setJintegralSIF(Jcp * crackData->getdomain_types() + i, 0.0);

            for (int j = 0; j < MODE; j++) 
            {
                setMintegralScalar(Jcp * (MODE * crackData->getdomain_types()) + i * MODE + j, 0.0);
                setMintegralSIF(Jcp * (MODE * crackData->getdomain_types()) + i * MODE + j, 0.0);
            }
        }

        // calculate the J-integral point parameters
        cout << endl << endl << "Jcp: " << Jcp << endl;
        double rangeMax = static_cast<double>(crackData->getnum_edge_crack_patch());                        // き裂前縁方向のパッチ数
        double crackJintegralNatural = (Jcp + 1.0) * (rangeMax / (crackData->getJintegralPoints() + 1.0));  // 自然座標系におけるη方向き裂先端座標[0 ~ き裂前縁パッチ数分の自然数]

        std::cout << "crackJintegralNatural [0, " << crackData->getnum_edge_crack_patch() << "] : " << crackJintegralNatural << std::endl;   

        vector<int> domainPatchID(crackData->getnum_singular_patch());                                                             // 計算パッチIDのセット
        vector<double> Jcenter(info->DIMENSION, 0.0);                                                                              // J 計算点η座標
        double JintegralArea = IIM_3D::makeCrackFrontVector(crackJintegralNatural, Jcp, domainPatchID, Jcenter, crackData, info);  // JintegralArea:積分領域の仮想き裂進展面積ΔA
        double degree = std::atan2(Jcenter[1], Jcenter[0]) * 180.0 / M_PI;                                                         // き裂前縁角度[°]
        std::cout << "crack leading edge angle: " << degree << std::endl;

        for (int ds = 0; ds < DOMAIN_SIDES; ds++) {                // 計算点のη方向両側の領域に対して処理
            double tempCenterNatural_eta = (ds == 0) ? crackJintegralNatural - crackData->getIntegralDomainNaturalWidth() : crackJintegralNatural;  // η方向積分領域の始点の設定

            for (int dp = 0; dp < crackData->getnum_singular_patch(); dp++) {        // 計算点を含む特異パッチ分の処理

                int patchID = domainPatchID[dp];  // 計算パッチIDの取得
                if (ds == 0)    cout << "processing patchID: " << patchID << endl;
                for (int divXi = 0; divXi < crackData->getIntegralDomainDivisions()[0]; divXi++) {                // ξ方向積分領域分割数
                    for (int divEta = 0; divEta < crackData->getIntegralDomainDivisions()[1]; divEta++) {         // η方向積分領域分割数
                        for (int divZeta = 0; divZeta < crackData->getIntegralDomainDivisions()[2]; divZeta++){   // ζ方向積分領域分割数   該当パッチの持つ要素番号によるループにしたい

                            vector<double> JTotalTerm1(crackData->getdomain_types() * crackData->getJintegralPoints(), 0.0), 
                                           JTotalTerm2(crackData->getdomain_types() * crackData->getJintegralPoints(), 0.0), 
                                           JTotalTerm3(crackData->getdomain_types() * crackData->getJintegralPoints(), 0.0);  //
                            vector<double> EMT1(MODE * crackData->getdomain_types() * crackData->getJintegralPoints(), 0.0), 
                                           EMT2(MODE * crackData->getdomain_types() * crackData->getJintegralPoints(), 0.0), 
                                           EMT3(MODE * crackData->getdomain_types() * crackData->getJintegralPoints(), 0.0);  //

                                        // クラスの定義
                                        PQSetting_3D actOverlay(info->DIMENSION);
                            
                                           // ガウス積分
                            for (int gaussID = 0; gaussID < GP_ON_ELEMENT; gaussID++) 
                            {
                                int integralElementID = 0;
                                double Jac = 0.0;
                                std::unique_ptr<QGradientFunction> qGrad = std::make_unique<QGradientFunction>(info->DIMENSION, *crackData);
                                vector<double> integralPointNatural(info->DIMENSION), integralPointNaturalTilde(info->DIMENSION), integralPointCoords(info->DIMENSION, 0.0);

                                //ガウス点をξ, η, ζ[-1, 1]座標から[0, 1]座標へ変換
                                integralPointNatural[0] = (static_cast<double>(divXi) / crackData->getIntegralDomainDivisions()[0] + (info->Gxi[gaussID * info->DIMENSION + 0] + 1.0) / (2.0 * crackData->getIntegralDomainDivisions()[0]));   
                                integralPointNatural[1] = tempCenterNatural_eta + (crackData->getIntegralDomainNaturalWidth() * (info->Gxi[gaussID * info->DIMENSION + 1] + 1.0)) / 2.0;
                                integralPointNatural[2] = (static_cast<double>(divZeta) / crackData->getIntegralDomainDivisions()[2] + (info->Gxi[gaussID * info->DIMENSION + 2] + 1.0) / (2.0 * crackData->getIntegralDomainDivisions()[2]));
                                
                                // ガウス点対するヤコビアン，q 関数の勾配を計算
                                IIM_3D::makeJintegralDomain(patchID, crackJintegralNatural, integralElementID, integralPointNatural, integralPointNaturalTilde, integralPointCoords, Jac, qGrad.get(), crackData, info); 

                                // calculate physical quantity in actual fields 
                                //std::unique_ptr<PQSetting_3D> actOverlay = std::make_unique<PQSetting_3D>(info->DIMENSION);
                                //computeActualField(gaussID, integralElementID, integralPointNaturalTilde, actOverlay.get(), info);
                                computeActualField(gaussID, integralElementID, integralPointNaturalTilde, &actOverlay, info);

                                const auto& StrainEnergyDensity = actOverlay.getStrainEnergyDensity();
                                const auto& Stress = actOverlay.getStress();
                                const auto& DispGrad1 = actOverlay.getDispGrad1();
                                const auto& DispGrad2 = actOverlay.getDispGrad2();
                                const auto& DispGrad3 = actOverlay.getDispGrad3();


                                // J積分
                                for (int domainSize = 0; domainSize < crackData->getdomain_types(); domainSize++)  
                                {
                                    vector<double> JfirstTerm1(info->DIMENSION, 0.0), JfirstTerm2(info->DIMENSION, 0.0), JfirstTerm3(info->DIMENSION, 0.0);     // J積分の第1項
                                    vector<double> JsecondTerm1(info->DIMENSION, 0.0), JsecondTerm2(info->DIMENSION, 0.0), JsecondTerm3(info->DIMENSION, 0.0);  // J積分の第2項

                                    for (int i = 0; i < info->DIMENSION; i++) 
                                    {
                                        JfirstTerm1[i] = StrainEnergyDensity * delta[0 * info->DIMENSION + i];
                                        JfirstTerm2[i] = StrainEnergyDensity * delta[1 * info->DIMENSION + i];
                                        JfirstTerm3[i] = StrainEnergyDensity * delta[2 * info->DIMENSION + i];
                                        // printf("firstTerm1[%d]: %.15e\n", i, JfirstTerm1[i]);
                                        // printf("firstTerm2[%d]: %.15e\n", i, JfirstTerm2[i]);
                                        // printf("firstTerm3[%d]: %.15e\n", i, JfirstTerm3[i]);   
                                    }

                                    // 第2項：応力・変位勾配項 (σᵢⱼ∂uⱼ/∂x₁)
                                    // x方向成分の行
                                    JsecondTerm1[0] =Stress[0] * DispGrad1[0] + 
                                                     Stress[3] * DispGrad1[1] + 
                                                     Stress[5] * DispGrad1[2];
                                    JsecondTerm1[1] =Stress[3] * DispGrad1[0] + 
                                                     Stress[1] * DispGrad1[1] + 
                                                     Stress[4] * DispGrad1[2];
                                    JsecondTerm1[2] =Stress[5] * DispGrad1[0] + 
                                                     Stress[4] * DispGrad1[1] + 
                                                     Stress[2] * DispGrad1[2];

                                    // y方向成分の行
                                    JsecondTerm2[0] =Stress[0] * DispGrad1[0] + 
                                                     Stress[3] * DispGrad1[1] + 
                                                     Stress[5] * DispGrad1[2];
                                    JsecondTerm2[1] =Stress[3] * DispGrad1[0] + 
                                                     Stress[1] * DispGrad1[1] + 
                                                     Stress[4] * DispGrad1[2];
                                    JsecondTerm2[2] =Stress[5] * DispGrad1[0] + 
                                                     Stress[4] * DispGrad1[1] + 
                                                     Stress[2] * DispGrad1[2];

                                    // z方向成分の行
                                    JsecondTerm3[0] =Stress[0] * DispGrad1[0] + 
                                                     Stress[3] * DispGrad1[1] + 
                                                     Stress[5] * DispGrad1[2];
                                    JsecondTerm3[1] =Stress[3] * DispGrad1[0] + 
                                                     Stress[1] * DispGrad1[1] + 
                                                     Stress[4] * DispGrad1[2];
                                    JsecondTerm3[2] =Stress[5] * DispGrad1[0] + 
                                                     Stress[4] * DispGrad1[1] + 
                                                     Stress[2] * DispGrad1[2];

                                    // for debug
                                    // for (int i = 0; i < info->DIMENSION; i++) 
                                    // {
                                    //     printf("domainSize : %d, secondTerm1[%d]: %.15e\n", domainSize, i, JsecondTerm1[i]);
                                    //     printf("domainSize : %d, secondTerm2[%d]: %.15e\n", domainSize, i, JsecondTerm2[i]);
                                    //     printf("domainSize : %d, secondTerm3[%d]: %.15e\n", domainSize, i, JsecondTerm3[i]);   
                                    // }

                                    #if 0 
                                    // J = -1/ΔA ∫_V (Wδ₁ᵢ - σᵢⱼ∂uⱼ/∂x₁)∂q/∂xᵢ dV`
                                    JTotalTerm1[domainSize] -= (JfirstTerm1[0] - JsecondTerm1[0]) * qGrad->getQGradient1()[domainSize * info->DIMENSION + 0] * info->w[gaussID] * Jac / JintegralArea;  
                                    JTotalTerm1[domainSize] -= (JfirstTerm1[1] - JsecondTerm1[1]) * qGrad->getQGradient2()[domainSize * info->DIMENSION + 0] * info->w[gaussID] * Jac / JintegralArea;    
                                    JTotalTerm1[domainSize] -= (JfirstTerm1[2] - JsecondTerm1[2]) * qGrad->getQGradient3()[domainSize * info->DIMENSION + 0] * info->w[gaussID] * Jac / JintegralArea;  
                                    JTotalTerm2[domainSize] -= (JfirstTerm2[0] - JsecondTerm2[0]) * qGrad->getQGradient1()[domainSize * info->DIMENSION + 1] * info->w[gaussID] * Jac / JintegralArea;  
                                    JTotalTerm2[domainSize] -= (JfirstTerm2[1] - JsecondTerm2[1]) * qGrad->getQGradient2()[domainSize * info->DIMENSION + 1] * info->w[gaussID] * Jac / JintegralArea;   
                                    JTotalTerm2[domainSize] -= (JfirstTerm2[2] - JsecondTerm2[2]) * qGrad->getQGradient3()[domainSize * info->DIMENSION + 1] * info->w[gaussID] * Jac / JintegralArea;  
                                    JTotalTerm3[domainSize] -= (JfirstTerm3[0] - JsecondTerm3[0]) * qGrad->getQGradient1()[domainSize * info->DIMENSION + 2] * info->w[gaussID] * Jac / JintegralArea;  
                                    JTotalTerm3[domainSize] -= (JfirstTerm3[1] - JsecondTerm3[1]) * qGrad->getQGradient2()[domainSize * info->DIMENSION + 2] * info->w[gaussID] * Jac / JintegralArea;   
                                    JTotalTerm3[domainSize] -= (JfirstTerm3[2] - JsecondTerm3[2]) * qGrad->getQGradient3()[domainSize * info->DIMENSION + 2] * info->w[gaussID] * Jac / JintegralArea;
                                    #endif
                                    
                                    int index = crackData->getdomain_types() * Jcp + domainSize;
                                    JTotalTerm1[index] -= (JfirstTerm1[0] - JsecondTerm1[0]) * qGrad->getQGradient1()[domainSize * info->DIMENSION + 0] * info->w[gaussID] * Jac / JintegralArea;  
                                    JTotalTerm1[index] -= (JfirstTerm1[1] - JsecondTerm1[1]) * qGrad->getQGradient2()[domainSize * info->DIMENSION + 0] * info->w[gaussID] * Jac / JintegralArea;    
                                    JTotalTerm1[index] -= (JfirstTerm1[2] - JsecondTerm1[2]) * qGrad->getQGradient3()[domainSize * info->DIMENSION + 0] * info->w[gaussID] * Jac / JintegralArea;  
                                    JTotalTerm2[index] -= (JfirstTerm2[0] - JsecondTerm2[0]) * qGrad->getQGradient1()[domainSize * info->DIMENSION + 1] * info->w[gaussID] * Jac / JintegralArea;  
                                    JTotalTerm2[index] -= (JfirstTerm2[1] - JsecondTerm2[1]) * qGrad->getQGradient2()[domainSize * info->DIMENSION + 1] * info->w[gaussID] * Jac / JintegralArea;   
                                    JTotalTerm2[index] -= (JfirstTerm2[2] - JsecondTerm2[2]) * qGrad->getQGradient3()[domainSize * info->DIMENSION + 1] * info->w[gaussID] * Jac / JintegralArea;  
                                    JTotalTerm3[index] -= (JfirstTerm3[0] - JsecondTerm3[0]) * qGrad->getQGradient1()[domainSize * info->DIMENSION + 2] * info->w[gaussID] * Jac / JintegralArea;  
                                    JTotalTerm3[index] -= (JfirstTerm3[1] - JsecondTerm3[1]) * qGrad->getQGradient2()[domainSize * info->DIMENSION + 2] * info->w[gaussID] * Jac / JintegralArea;   
                                    JTotalTerm3[index] -= (JfirstTerm3[2] - JsecondTerm3[2]) * qGrad->getQGradient3()[domainSize * info->DIMENSION + 2] * info->w[gaussID] * Jac / JintegralArea;

                                    // for debug
                                    // for (int i = 0; i < info->DIMENSION; i++) 
                                    // {
                                    //     printf("domainSize : %d, JTotalTerm1[%d]: %.15e\n", domainSize, i, JTotalTerm1[i]);
                                    //     printf("domainSize : %d, JTotalTerm2[%d]: %.15e\n", domainSize, i, JTotalTerm2[i]);
                                    //     printf("domainSize : %d, JTotalTerm3[%d]: %.15e\n", domainSize, i, JTotalTerm3[i]);
                                    // }

                                    // int index = domainSize;
                                    setJintegralScalar(index, getJintegralScalar()[index] + 
                                    JTotalTerm1[index] * crackFrontVector[2 * info->DIMENSION + 0] + 
                                    JTotalTerm2[index] * crackFrontVector[2 * info->DIMENSION + 1] + 
                                    JTotalTerm3[index] * crackFrontVector[2 * info->DIMENSION + 2]);
                                    // for debug
                                    // printf("JintegralScalar[%d]: %.15e\n", domainSize, getJintegralScalar()[domainSize]);
                                }

                                // IIM

                                // auxの定義
                                PQSetting_3D aux(info->DIMENSION);

                                for (int mode = 0; mode < MODE; mode++) // 計算する破壊モードの数だけ処理
                                {
                                    // calculate physical quantity auxiliary fields
                                    IIM_3D::computeAuxiliaryField(mode + 1, Jcenter, integralPointCoords, auxField.get(), &actOverlay, &aux, info);

                                const auto& StrainEnergyDensity_aux = aux.getStrainEnergyDensity();
                                const auto& Stress_aux = aux.getStress();
                                const auto& DispGrad1_aux = aux.getDispGrad1();
                                const auto& DispGrad2_aux = aux.getDispGrad2();
                                const auto& DispGrad3_aux = aux.getDispGrad3();


                                    // EMT = ∫_V W(1,2)δ_1i - (σ_ij(1) * ∂u_j(2) / ∂x_1 + σ_ij(2) * ∂u_j(1) / ∂x_1 ) ∂q/∂x_i dV
                                    for (int domainSize = 0; domainSize < crackData->getdomain_types(); domainSize++)  // 
                                    {
                                        vector<double> firstEMT1(info->DIMENSION), secondEMT1(info->DIMENSION), thirdEMT1(info->DIMENSION);  // EMT x
                                        vector<double> firstEMT2(info->DIMENSION), secondEMT2(info->DIMENSION), thirdEMT2(info->DIMENSION);  // EMT y
                                        vector<double> firstEMT3(info->DIMENSION), secondEMT3(info->DIMENSION), thirdEMT3(info->DIMENSION);  // EMT z

                                        // first term (strain energy density * delta)
                                        for (int i = 0; i < info->DIMENSION; i++) 
                                        {
                                            firstEMT1[i] = StrainEnergyDensity_aux * delta[0 * info->DIMENSION + i];  // EMT 第一項 x
                                            firstEMT2[i] = StrainEnergyDensity_aux * delta[1 * info->DIMENSION + i];  // EMT 第一項 y
                                            firstEMT3[i] = StrainEnergyDensity_aux * delta[2 * info->DIMENSION + i];  // EMT 第一項 z
                                        }

                                        // second term (actual stress * auxiliary displacement gradient)
                                        // x方向成分
                                        secondEMT1[0] =Stress[0] * DispGrad1_aux[0] +Stress[3] * DispGrad1_aux[1] +Stress[5] * DispGrad1_aux[2];  // σxx * ∂ux/∂x + σxy * ∂ux/∂y + σzx * ∂ux/∂z (j = 1)
                                        secondEMT1[1] =Stress[3] * DispGrad1_aux[0] +Stress[1] * DispGrad1_aux[1] +Stress[4] * DispGrad1_aux[2];  // σxy * ∂ux/∂x + σyy * ∂ux/∂y + σyz * ∂ux/∂z (j = 2)
                                        secondEMT1[2] =Stress[5] * DispGrad1_aux[0] +Stress[4] * DispGrad1_aux[1] +Stress[2] * DispGrad1_aux[2];  // σxz * ∂ux/∂x + σyz * ∂ux/∂y + σzz * ∂ux/∂z (j = 3)

                                        // y方向成分
                                        secondEMT2[0] =Stress[0] * DispGrad2_aux[0] +Stress[3] * DispGrad2_aux[1] +Stress[5] * DispGrad2_aux[2];  // σxx * ∂uy/∂x + σxy * ∂uy/∂y + σxz * ∂uy/∂z (j = 1)
                                        secondEMT2[1] =Stress[3] * DispGrad2_aux[0] +Stress[1] * DispGrad2_aux[1] +Stress[4] * DispGrad2_aux[2];  // σxy * ∂uy/∂x + σyy * ∂uy/∂y + σyz * ∂uy/∂z (j = 2)
                                        secondEMT2[2] =Stress[5] * DispGrad2_aux[0] +Stress[4] * DispGrad2_aux[1] +Stress[2] * DispGrad2_aux[2];  // σxz * ∂uy/∂x + σyz * ∂uy/∂y + σzz * ∂uy/∂z (j = 3)

                                        // z方向成分
                                        secondEMT3[0] =Stress[0] * DispGrad3_aux[0] +Stress[3] * DispGrad3_aux[1] +Stress[5] * DispGrad3_aux[2];  // σxx * ∂uz/∂x + σxy * ∂uz/∂y + σxz * ∂uz/∂z (j = 1)
                                        secondEMT3[1] =Stress[3] * DispGrad3_aux[0] +Stress[1] * DispGrad3_aux[1] +Stress[4] * DispGrad3_aux[2];  // σxy * ∂uz/∂x + σyy * ∂uz/∂y + σyz * ∂uz/∂z (j = 2)
                                        secondEMT3[2] =Stress[5] * DispGrad3_aux[0] +Stress[4] * DispGrad3_aux[1] +Stress[2] * DispGrad3_aux[2];  // σxz * ∂uz/∂x + σyz * ∂uz/∂y + σzz * ∂uz/∂z (j = 3)

                                        // third term (auxiliary stress * actual displacement gradient)
                                        thirdEMT1[0] = Stress_aux[0] * DispGrad1[0] + Stress_aux[3] * DispGrad1[1] + Stress_aux[5] * DispGrad1[2];  // σxx * ∂ux/∂x + σxy * ∂ux/∂y + σxz * ∂ux/∂z (j = 1)
                                        thirdEMT1[1] = Stress_aux[3] * DispGrad1[0] + Stress_aux[1] * DispGrad1[1] + Stress_aux[4] * DispGrad1[2];  // σxy * ∂ux/∂x + σyy * ∂ux/∂y + σyz * ∂ux/∂z (j = 2)
                                        thirdEMT1[2] = Stress_aux[5] * DispGrad1[0] + Stress_aux[4] * DispGrad1[1] + Stress_aux[2] * DispGrad1[2];  // σxz * ∂ux/∂x + σyz * ∂ux/∂y + σzz * ∂ux/∂z (j = 3)

                                        thirdEMT2[0] = Stress_aux[0] * DispGrad2[0] + Stress_aux[3] * DispGrad2[1] + Stress_aux[5] * DispGrad2[2];  // σxx * ∂uy/∂x + σxy * ∂uy/∂y + σxz * ∂uy/∂z (j = 1)
                                        thirdEMT2[1] = Stress_aux[3] * DispGrad2[0] + Stress_aux[1] * DispGrad2[1] + Stress_aux[4] * DispGrad2[2];  // σxy * ∂uy/∂x + σyy * ∂uy/∂y + σyz * ∂uy/∂z (j = 2)
                                        thirdEMT2[2] = Stress_aux[5] * DispGrad2[0] + Stress_aux[4] * DispGrad2[1] + Stress_aux[2] * DispGrad2[2];  // σxz * ∂uy/∂x + σyz * ∂uy/∂y + σzz * ∂uy/∂z (j = 3)

                                        thirdEMT3[0] = Stress_aux[0] * DispGrad3[0] + Stress_aux[3] * DispGrad3[1] + Stress_aux[5] * DispGrad3[2];  // σxx * ∂uz/∂x + σxy * ∂uz/∂y + σxz * ∂uz/∂z (j = 1)
                                        thirdEMT3[1] = Stress_aux[3] * DispGrad3[0] + Stress_aux[1] * DispGrad3[1] + Stress_aux[4] * DispGrad3[2];  // σxy * ∂uz/∂x + σyy * ∂uz/∂y + σyz * ∂uz/∂z (j = 2)
                                        thirdEMT3[2] = Stress_aux[5] * DispGrad3[0] + Stress_aux[4] * DispGrad3[1] + Stress_aux[2] * DispGrad3[2];  // σxz * ∂uz/∂x + σyz * ∂uz/∂y + σzz * ∂uz/∂z (j = 3)

                                        #if 0
                                        // total EMT
                                        EMT1[mode * DOMAIN_TYPES + domainSize] -= (firstEMT1[0] - secondEMT1[0] - thirdEMT1[0]) * qGrad->getQGradient1()[domainSize * info->DIMENSION + 0] * info->w[gaussID] * Jac / JintegralArea;
                                        EMT1[mode * DOMAIN_TYPES + domainSize] -= (firstEMT1[1] - secondEMT1[1] - thirdEMT1[1]) * qGrad->getQGradient2()[domainSize * info->DIMENSION + 0] * info->w[gaussID] * Jac / JintegralArea;
                                        EMT1[mode * DOMAIN_TYPES + domainSize] -= (firstEMT1[2] - secondEMT1[2] - thirdEMT1[2]) * qGrad->getQGradient3()[domainSize * info->DIMENSION + 0] * info->w[gaussID] * Jac / JintegralArea;

                                        EMT2[mode * DOMAIN_TYPES + domainSize] -= (firstEMT2[0] - secondEMT2[0] - thirdEMT2[0]) * qGrad->getQGradient1()[domainSize * info->DIMENSION + 1] * info->w[gaussID] * Jac / JintegralArea;
                                        EMT2[mode * DOMAIN_TYPES + domainSize] -= (firstEMT2[1] - secondEMT2[1] - thirdEMT2[1]) * qGrad->getQGradient2()[domainSize * info->DIMENSION + 1] * info->w[gaussID] * Jac / JintegralArea;
                                        EMT2[mode * DOMAIN_TYPES + domainSize] -= (firstEMT2[2] - secondEMT2[2] - thirdEMT2[2]) * qGrad->getQGradient3()[domainSize * info->DIMENSION + 1] * info->w[gaussID] * Jac / JintegralArea;

                                        EMT3[mode * DOMAIN_TYPES + domainSize] -= (firstEMT3[0] - secondEMT3[0] - thirdEMT3[0]) * qGrad->getQGradient1()[domainSize * info->DIMENSION + 2] * info->w[gaussID] * Jac / JintegralArea;
                                        EMT3[mode * DOMAIN_TYPES + domainSize] -= (firstEMT3[1] - secondEMT3[1] - thirdEMT3[1]) * qGrad->getQGradient2()[domainSize * info->DIMENSION + 2] * info->w[gaussID] * Jac / JintegralArea;
                                        EMT3[mode * DOMAIN_TYPES + domainSize] -= (firstEMT3[2] - secondEMT3[2] - thirdEMT3[2]) * qGrad->getQGradient3()[domainSize * info->DIMENSION + 2] * info->w[gaussID] * Jac / JintegralArea;
                                        #endif 
                                        int index = Jcp * MODE * crackData->getdomain_types() + mode * crackData->getdomain_types() + domainSize;
                                        EMT1[index] -= (firstEMT1[0] - secondEMT1[0] - thirdEMT1[0]) * qGrad->getQGradient1()[domainSize * info->DIMENSION + 0] * info->w[gaussID] * Jac / JintegralArea;
                                        EMT1[index] -= (firstEMT1[1] - secondEMT1[1] - thirdEMT1[1]) * qGrad->getQGradient2()[domainSize * info->DIMENSION + 0] * info->w[gaussID] * Jac / JintegralArea;
                                        EMT1[index] -= (firstEMT1[2] - secondEMT1[2] - thirdEMT1[2]) * qGrad->getQGradient3()[domainSize * info->DIMENSION + 0] * info->w[gaussID] * Jac / JintegralArea;

                                        EMT2[index] -= (firstEMT2[0] - secondEMT2[0] - thirdEMT2[0]) * qGrad->getQGradient1()[domainSize * info->DIMENSION + 1] * info->w[gaussID] * Jac / JintegralArea;
                                        EMT2[index] -= (firstEMT2[1] - secondEMT2[1] - thirdEMT2[1]) * qGrad->getQGradient2()[domainSize * info->DIMENSION + 1] * info->w[gaussID] * Jac / JintegralArea;
                                        EMT2[index] -= (firstEMT2[2] - secondEMT2[2] - thirdEMT2[2]) * qGrad->getQGradient3()[domainSize * info->DIMENSION + 1] * info->w[gaussID] * Jac / JintegralArea;

                                        EMT3[index] -= (firstEMT3[0] - secondEMT3[0] - thirdEMT3[0]) * qGrad->getQGradient1()[domainSize * info->DIMENSION + 2] * info->w[gaussID] * Jac / JintegralArea;
                                        EMT3[index] -= (firstEMT3[1] - secondEMT3[1] - thirdEMT3[1]) * qGrad->getQGradient2()[domainSize * info->DIMENSION + 2] * info->w[gaussID] * Jac / JintegralArea;
                                        EMT3[index] -= (firstEMT3[2] - secondEMT3[2] - thirdEMT3[2]) * qGrad->getQGradient3()[domainSize * info->DIMENSION + 2] * info->w[gaussID] * Jac / JintegralArea;
                                    }
                                }
                            }
                            // き裂前縁座標への射影
                            for (int mode = 0; mode < MODE; mode++) 
                            {
                                for (int domainSize = 0; domainSize < crackData->getdomain_types(); domainSize++) 
                                {
                                    // int index = mode * DOMAIN_TYPES + domainSize;
                                    int index = Jcp * MODE * crackData->getdomain_types() + mode * crackData->getdomain_types() + domainSize;
                                    setMintegralScalar(index, getMintegralScalar()[index] + 
                                                       EMT1[index] * crackFrontVector[2 * info->DIMENSION + 0] + 
                                                       EMT2[index] * crackFrontVector[2 * info->DIMENSION + 1] + 
                                                       EMT3[index] * crackFrontVector[2 * info->DIMENSION + 2]);
                                }
                            }
                        }
                    }
                }
            }
        }

        //  モデルサイズによる積分値の調整
        double model_size_para = 0;
        // if (MODEL_SIZE_PARA_3D == 0) 
        if (crackData->getmodel_size_para_3D() == 0) 
        {   
            printf("\nModel both sides of the crack surface : \n");
            model_size_para = 2.0;
        }
        if (crackData->getmodel_size_para_3D() == 1 || crackData->getmodel_size_para_3D() == 2)
        {   
            printf("\nAnalysis Model : Model with only one side of the crack surface\n");
            model_size_para = 1.0;
        }

        // Jintegral to SIF
        for (int domainSize = 0; domainSize < crackData->getdomain_types(); domainSize++) 
        {   
            int index = crackData->getdomain_types() * Jcp + domainSize;
            setJintegralSIF(domainSize, sqrt(getJintegralScalar()[domainSize] * E / (1.0 - nu * nu) * model_size_para));  // 平面ひずみ条件における変換式
            setJintegralSIF(index, sqrt(getJintegralScalar()[index] * E / (1.0 - nu * nu) * model_size_para));  // 平面ひずみ条件における変換式
            printf("\n******************************************************************************\n");
            cout << "size: all_element  - " << domainSize << endl;
            printf("Jintegral = %.15e  SIFs = %.15e\n", getJintegralScalar()[domainSize], getJintegralSIF()[domainSize]);
            printf("******************************************************************************\n");
        }

        // Mintegral to SIFs
        for (int mode = 0; mode < MODE; mode++) 
        {
            for (int domainSize = 0; domainSize < crackData->getdomain_types(); domainSize++) {
                // int index = mode * DOMAIN_TYPES + domainSize;
                int index = Jcp * MODE * crackData->getdomain_types() + mode * crackData->getdomain_types() + domainSize;
                if (mode != 2) {  // mode I, II
                    setMintegralSIF(index, 
                        getMintegralScalar()[index] * E / (model_size_para * (1.0 - nu * nu)));
                } else {         // mode III
                    setMintegralSIF(index, 
                        getMintegralScalar()[index] * E / (model_size_para * (1.0 + nu)));
                }

                printf("\n******************************************************************************\n");
                cout << "mode: " << mode + 1 << ", size: all_element  - " << domainSize << endl;
                printf("Mintegral = %.15e  SIF = %.15e\n", getMintegralScalar()[index], getMintegralSIF()[index]);
                printf("******************************************************************************\n");
            }
        }
        IIM_3D::writeSIFsToCSV("SIFs.csv", Jcp, Jcenter, degree, info, crackData);
    }
}


// き裂前縁に関する情報の計算
// 計算点のパッチIDの取得，各空間における座標，積分領域におけるき裂長さ仮想き裂進展面積の計算
double IIM_3D::makeCrackFrontVector(double &crackJintegralNatural, int Jcp, vector<int> &domainPatchID, vector<double> &Jcenter, CrackDataSetting_3D *crackData, information *info)
{
    int i, j, k,divideID;
    int local_init_patch = info->Total_Patch_to_mesh[Total_mesh - 1];
    // int local_init_patch = info->Total_Patch_to_mesh[1];  

    // き裂前縁方向のパッチインデックスの取得
    int flag = 0;        
    for (i = 0; i < crackData->getnum_edge_crack_patch(); i++) 
    {   
        if (crackJintegralNatural <= static_cast<double>(i + 1)) 
        {   
            flag = i;
            crackJintegralNatural = crackJintegralNatural - static_cast<double>(i);
            break;
        }
    }
    
    // 計算パッチのセットを作成
    for (j = 0; j < crackData->getnum_singular_patch(); j++) 
    {
        domainPatchID[j] = crackData->getpatchList()[(crackData->getnum_singular_patch() * flag + j)] + local_init_patch;
    }
    // for debug
    // for (j = 0; j < crackData->getnum_singular_patch(); j++){
    //     cout << "domainPatchID[" << j << "] : " << domainPatchID[j] << endl;
    // }

    double tempEdge_start = crackJintegralNatural - crackData->getIntegralDomainNaturalWidth();  // η 方向計算点積分領域始点 = η[0, 1]座標における計算点 - η 方向積分領域幅
    double tempEdge_end = crackJintegralNatural + crackData->getIntegralDomainNaturalWidth();    // η 方向計算点積分領域終点 = η[0, 1]座標における計算点 + η 方向積分領域幅

    // η方向計算点積分領域がパッチをまたぐ場合の処理
    if (MOVE_DOMAIN == 1)
    {
        const double epsilon = 1e-12;
        const int maxIteration = 100;
        int iteration = 0;

        while (tempEdge_start < 0.0 - epsilon || tempEdge_end > 1.0 + epsilon)
        {   
            // η方向計算点積分領域始点側がパッチをまたぐ場合の処理
            if (tempEdge_start < 0.0 - epsilon)  
            {
                cout << "Jcp : " << Jcp << " ---Integral area spans multiple patches" << endl;
                cout << "Moves the calculation point in the positive direction of the η-coordinate" << endl;
                crackJintegralNatural += crackData->getIntegralDomainNaturalWidth() / 10;
            }

            // η方向計算点積分領域終点側がパッチをまたぐ場合の処理
            if (tempEdge_end > 1.0 + epsilon) 
            {
                cout << "Jcp : " << Jcp << " ---Integral area spans multiple patches" << endl;
                cout << "Moves the calculation point in the negative direction of the η-coordinate" << endl;
                crackJintegralNatural -= crackData->getIntegralDomainNaturalWidth() / 10;
            }

            // 新しい積分領域の始点と終点
            tempEdge_start = crackJintegralNatural - crackData->getIntegralDomainNaturalWidth();
            tempEdge_end = crackJintegralNatural + crackData->getIntegralDomainNaturalWidth();

            iteration++;
            if (iteration > maxIteration) 
            {
                std::cerr << "Error: The number of iterations exceeded the maximum number of iterations." << std::endl;
                exit(1);
            }
        }
    }
    crackEdgeJcpNatural.push_back(crackJintegralNatural);
    // std::cout << "crackJintegralNatural [0, 1]: " << crackJintegralNatural << std::endl;
    // std::cout << "crackEdgeNatural: " << getcrackEdgeJcpNatural()[Jcp] << std::endl;

    // find the center location
    vector<double> tempCenterNatural(info->DIMENSION, 0.0), tempCenterNaturalTilde(info->DIMENSION, 0.0);
    int tempElementID = 0;
    tempCenterNatural[1] = crackJintegralNatural;  // [0.0, crackJintegralNatural, 0.0]

    for (i = 0; i < elements_in_patch[domainPatchID[0]]; i++)  // 各パッチに含まれる要素数分の処理
    {  
        if (IIM_3D::isCoordinateInElement(tempCenterNatural, domainPatchID[0], elementID_in_patch[domainPatchID[0]][i], info))  // tempCenterNaturalがelementID_in_patch[domainPatchID[0]][i]どの要素内にあるか判定して，該当要素のみ処理
        {
            tempElementID = elementID_in_patch[domainPatchID[0]][i];                                                                     // 該当要素IDを格納
            tempCenterNaturalTilde = math_->convertNaturalCoordinatesToTilde(domainPatchID[0], tempElementID, tempCenterNatural, info);  // 計算点を[0, 1]座標から[-1, 1]座標に変換
            break;
        }
    }
    crackEdgeJcpElementID.push_back(tempElementID);
    // std::cout << "tempElementID in tempNatural: " << tempElementID << std::endl;
    // std::cout << "crackEdgeJcpElementID: " << getcrackEdgeJcpElementID()[Jcp] << std::endl;

    // calculate the center location in physical coordinates
    // Jcenter = Σ Ni * xi
    for (i = 0; i < info->No_Control_point_ON_ELEMENT[domainPatchID[0]]; i++) {
        for (j = 0; j < info->DIMENSION; j++) {
            Jcenter[j] += Shape_func(i, tempCenterNaturalTilde.data(), tempElementID, info) * info->Node_Coordinate[info->Controlpoint_of_Element[tempElementID * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + j];
        }
    } 

    for (i = 0; i < info->DIMENSION; i++) 
    {   
        int index = Jcp * info->DIMENSION + i;
        setCalculationPointCoords(index, Jcenter[i]);
        // std::cout << "Jcenter in physical coordinate[" << i << "]: " << Jcenter[i] << std::endl;
        std::cout << "Jcenter in physical coordinate[" << i << "]: " << getCalculationPointCoords()[index] << std::endl;
    }


    // calculate virtual crack extension value
    double length = 0.0;         // 計算領域のき裂長さ
    double JintegralArea = 0.0;  // 計算領域の仮想き裂進展面積
    vector<double> tempEdgeNatural(info->DIMENSION, 0.0);                    
    vector<double> tempEdgeNaturalTilde(info->DIMENSION, 0.0);               
    vector<double> pointCoords(crackData->getdomain_area_divisions() * info->DIMENSION, 0.0);  
    vector<double> G1d(crackData->getdomain_area_divisions(), 0.0);

    for (divideID = 0; divideID < crackData->getdomain_area_divisions(); divideID++)  // η方向の積分領域の分割数だけ処理
    {
        G1d[divideID] = -1.0 + 2.0 / (crackData->getdomain_area_divisions() - 1.0) * divideID;  // η方向分割領域を[-1, 1]に変換
        // tempEdgeNatural = [0, 1]座標における各分割領域に対する始点
        tempEdgeNatural[0] = 0.0;
        tempEdgeNatural[1] = tempEdge_start + (crackData->getIntegralDomainNaturalWidth() * (G1d[divideID] + 1.0));  // (G1d[divideID] + 1.0) : η方向分割領域を[0, 1]にとして扱うため
        tempEdgeNatural[2] = 0.0;

        for (i = 0; i < elements_in_patch[domainPatchID[0]]; i++) // 各パッチの要素数分の処理
        {
            if (IIM_3D::isCoordinateInElement(tempEdgeNatural, domainPatchID[0], elementID_in_patch[domainPatchID[0]][i], info))  // tempEdgeNaturalがelementID_in_patch[patchID][i]どの要素内にあるか判定して，該当要素のみ処理
            {
                tempElementID = elementID_in_patch[domainPatchID[0]][i];  // 該当要素IDを抽出
                tempEdgeNaturalTilde = math_->convertNaturalCoordinatesToTilde(domainPatchID[0], tempElementID, tempEdgeNatural, info);  // tempEdgeNaturalを[0, 1]座標から[-1, 1]座標に変換

                // calculate the J-domian edge location in physical coordinates
                // pointCoords = Σ Ni * xi
                for (j = 0; j < info->No_Control_point_ON_ELEMENT[domainPatchID[0]]; j++) 
                {
                    for (k = 0; k < info->DIMENSION; k++) {
                        pointCoords[divideID * info->DIMENSION + k] += Shape_func(j, tempEdgeNaturalTilde.data(), tempElementID, info) * info->Node_Coordinate[info->Controlpoint_of_Element[tempElementID * MAX_NO_CP_ON_ELEMENT + j] * (info->DIMENSION + 1) + k];
                    }
                }

                // 計算領域のき裂長さを計算
                if (divideID > 0) length += sqrt((pointCoords[(divideID - 1) * info->DIMENSION + 0] - pointCoords[divideID * info->DIMENSION + 0]) * (pointCoords[(divideID - 1) * info->DIMENSION + 0] - pointCoords[divideID * info->DIMENSION + 0]) +
                                                 (pointCoords[(divideID - 1) * info->DIMENSION + 1] - pointCoords[divideID * info->DIMENSION + 1]) * (pointCoords[(divideID - 1) * info->DIMENSION + 1] - pointCoords[divideID * info->DIMENSION + 1]) +
                                                 (pointCoords[(divideID - 1) * info->DIMENSION + 2] - pointCoords[divideID * info->DIMENSION + 2]) * (pointCoords[(divideID - 1) * info->DIMENSION + 2] - pointCoords[divideID * info->DIMENSION + 2]));
                break;
            }
        }
    }
    JintegralArea = length / 2.0;  // 仮想き裂面積 ΔA = length * q / 2 (q = 1)
    
    // き裂長さと面積のエラーハンドリング
    if (length < 0 || JintegralArea < 0) 
    {
        std::cerr << "Error: The crack length or the virtual crack extension area is negative." << std::endl;
        exit(1);
    }
    std::cout << "Crack Length: " << length << std::endl;
    std::cout << "Jintegral Area: " << JintegralArea << std::endl;


    // set crack front vector
    vector<double> tempVector1(info->DIMENSION), tempVector2(info->DIMENSION);
    for (i = 0; i < info->DIMENSION; i++) 
    {
        tempVector1[i] = pointCoords[0 * info->DIMENSION + i] - Jcenter[i];                          // 計算領域始点からき裂中心(J 計算点)までのベクトル
        tempVector2[i] = Jcenter[i] - pointCoords[(crackData->getdomain_area_divisions() - 1) * info->DIMENSION + i];  // き裂中心(J 計算点)から計算領域終点までのベクトル
    }


    IIM_3D::makeOrthogonalBasisVector(tempVector1, tempVector2, info);  // き裂前縁法線，接線およびき裂面法線の基底ベクトルを計算
    IIM_3D::storeJcpCrackFrontVector(info, Jcp);                        // き裂前縁法線（き裂前縁ベクトル）を格納

    return JintegralArea;
}


void IIM_3D::readCrackFrontData(CrackDataSetting_3D *crackData, information *info)
{
    // iterate over local patch elements
    for (int re = info->real_Total_Element_to_mesh[Total_mesh - 1]; re < info->real_Total_Element_to_mesh[Total_mesh]; re++)  // 現在参照パッチの要素数分の処理
    {
        int e = info->real_element[re];
        elementID_in_patch[info->Element_patch[e]].push_back(e);  // 各パッチに属する要素ID のリストを追加  
        ++elements_in_patch[info->Element_patch[e]];              // 各パッチに含まれる要素の数をカウント   
    }

    for (int Jcp = 0; Jcp < crackData->getJintegralPoints(); Jcp++)  // き裂前縁方向に設定した計算点数分だけ処理
    {
        // calculate the J-integral point parameters
        cout << endl << endl << "Jcp: " << Jcp << endl;
        double rangeMax = static_cast<double>(crackData->getnum_edge_crack_patch());                                                         // き裂前縁方向のパッチ数
        double crackJintegralNatural = (Jcp + 1.0) * (rangeMax / (crackData->getJintegralPoints() + 1.0));                         // 自然座標系におけるη方向き裂先端座標[0 ~ き裂前縁パッチ数分の自然数]

        std::cout << "crackJintegralNatural [0, " << crackData->getnum_edge_crack_patch() << "] : " << crackJintegralNatural << std::endl;   

        vector<int> domainPatchID(crackData->getnum_singular_patch());                                                             // 計算パッチIDのセット
        vector<double> Jcenter(info->DIMENSION, 0.0);                                                                              // J 計算点η座標
        // double JintegralArea = IIM_3D::makeCrackFrontVector(crackJintegralNatural, Jcp, domainPatchID, Jcenter, crackData, info);  // JintegralArea:積分領域の仮想き裂進展面積ΔA
        IIM_3D::makeCrackFrontVector(crackJintegralNatural, Jcp, domainPatchID, Jcenter, crackData, info);  // JintegralArea:積分領域の仮想き裂進展面積ΔA
        double degree = std::atan2(Jcenter[1], Jcenter[0]) * 180.0 / M_PI;                                                         // き裂前縁角度[°]
        std::cout << "crack leading edge angle: " << degree << std::endl;
    }
}


void IIM_3D::storeJcpCrackFrontVector(information *info, int Jcp)
{
    // cout << "Jcp: " << Jcp << endl;
    for (int i = 0; i < info->DIMENSION; i++)  // 3種類のベクトル（binormal, tangent, normal）
    {   
        for (int dim = 0; dim < info->DIMENSION; dim++) // 各ベクトルの成分（x, y, z）
        {
            int index_output = Jcp * info->DIMENSION * info->DIMENSION + i * info->DIMENSION + dim;
            int index_input = i * info->DIMENSION + dim;
            setJcpCrackFrontVector(index_output, crackFrontVector[index_input]);
            // cout << "Crack Front Vector[" << Jcp << "][" << i << "][" << dim << "]: " << getJcpCrackFrontVector()[index_output] << endl;
            // cout << "Crack Front Vector[" << Jcp << "][" << i << "][" << dim << "]: " << crackFrontVector[index_output] << endl;
        }
    }
}


// き裂先端を含むパッチで、η方向にき裂先端を定義することを前提とした関数
void IIM_3D::makeJintegralDomain(const int patchID, const double crackJintegralNatural, int &integralElementID, vector<double> &integralPointNatural, vector<double> &integralPointNaturalTilde, vector<double> &integralPointCoords, double &integralJac, QGradientFunction *qGrad, CrackDataSetting_3D *crackData, information *info)
{
    int i, j, k, l;

    for (i = 0; i < elements_in_patch[patchID]; i++) {
        if (IIM_3D::isCoordinateInElement(integralPointNatural, patchID, elementID_in_patch[patchID][i], info))  // integralPointNatural(ガウス点)が elementID_in_patch[patchID][i] どの要素内にあるか判定して，該当要素のみ処理
        {  
            // ガウス点を自然座標[0, 1]から[-1, 1]へ変換
            integralElementID = elementID_in_patch[patchID][i];                                                                           // 該当要素IDを格納
            integralPointNaturalTilde = math_->convertNaturalCoordinatesToTilde(patchID, integralElementID, integralPointNatural, info);  // integralPointNaturalを[0, 1]座標から[-1, 1]座標に変換

            // ガウス点を[-1, 1]座標から物理座標へ変換, x = Nixi
            for (j = 0; j < info->No_Control_point_ON_ELEMENT[patchID]; j++) 
            {
                for (k = 0; k < info->DIMENSION; k++) {
                    integralPointCoords[k] += Shape_func(j, integralPointNaturalTilde.data(), integralElementID, info) * info->Node_Coordinate[info->Controlpoint_of_Element[integralElementID * MAX_NO_CP_ON_ELEMENT + j] * (info->DIMENSION + 1) + k]; 
                }
            }

            // differential coefficient (dXiTilde_dXiPrime : ∂ξ^~ / ∂ξ', dXiPrime_dXiPrimeTilde : ∂ξ' / ∂ξ'^~) ξ'^~:[-1, 1], ξ^~:[0, 1], ξ':[0, 1/N]  
            vector<double> dXiTilde_dXiPrime(info->DIMENSION), dXiPrime_dXiPrimeTilde(info->DIMENSION);
            // dξ^~ / dξ' = 2 / (ξ_{i+1} - ξ_i), η,ζ方向も同様に計算
            for (j = 0; j < info->DIMENSION; j++) 
            {
                dXiTilde_dXiPrime[j] = 2.0 / (info->Position_Knots[info->Total_Knot_to_patch_dim[patchID * info->DIMENSION + j] + info->Order[patchID * info->DIMENSION + j] + info->ENC[integralElementID * info->DIMENSION + j] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[patchID * info->DIMENSION + j] + info->Order[patchID * info->DIMENSION + j] + info->ENC[integralElementID * info->DIMENSION + j] + 0]);
            }

            dXiPrime_dXiPrimeTilde[0] = 1.0 / (2.0 * crackData->getIntegralDomainDivisions()[0]);  // ∂ξ' / ∂ξ`^~ = 1 / (2 * ξ 方向要素数)
            dXiPrime_dXiPrimeTilde[2] = 1.0 / (2.0 * crackData->getIntegralDomainDivisions()[2]);  // ∂ζ' / ∂ζ`^~ = 1 / (2 * ζ 方向要素数) 
            dXiPrime_dXiPrimeTilde[1] = crackData->getIntegralDomainNaturalWidth() / 2.0;          // ∂η' / ∂η`^~ = η方向積分領域幅 / 2


            // [∂x/∂ξ, ∂x/∂η, ∂x/∂ζ, ∂y/∂ξ, ∂y/∂η, ∂y/∂ζ, ∂z/∂ξ, ∂z/∂η, ∂z/∂ζ] を計算
            vector<double> a(info->DIMENSION * info->DIMENSION, 0.0);
            for (j = 0; j < info->DIMENSION; j++)
            {
                for (k = 0; k < info->DIMENSION; k++) 
                {
                    for (l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[integralElementID]]; l++) 
                    {
                        a[j * info->DIMENSION + k] += dShape_func(l, k, integralPointNaturalTilde.data(), integralElementID, info) * info->Node_Coordinate[info->Controlpoint_of_Element[integralElementID * MAX_NO_CP_ON_ELEMENT + l] * (info->DIMENSION + 1) + j];  // a[][] = ∂x/∂ξ
                    }
                }
            }

            // ヤコビアン計算
            vector<double> Jac(info->DIMENSION * info->DIMENSION);
            #pragma omp parallel for collapse(2) if(OpenMP)
            for(j = 0; j < info->DIMENSION; j++) 
            {
                for (k = 0; k < info->DIMENSION; k++) 
                {
                    Jac[j * info->DIMENSION + k] = a[j * info->DIMENSION + k] * dXiTilde_dXiPrime[k] * dXiPrime_dXiPrimeTilde[k];  //  Jac = (∂x/∂ξ) * (∂ξ^~ / ∂ξ') *  (∂ξ' / ∂ξ`^~) = ∂x / ∂ξ^~
                }
            }

            integralJac = InverseMatrix_3x3(Jac.data());  // ヤコビアンの逆行列 Jac^-1 = ∂ξ^~ / ∂x
            // calculate the q gradient function
            qGrad->makeQGradient(patchID, crackJintegralNatural, a, dXiTilde_dXiPrime, integralPointNatural, crackData, info);

            break;
        }
    }
}

void IIM_3D::computeActualField(const int gaussID, const int integralElementID, vector<double> &integralPointNaturalTilde, PQSetting_3D *actOverlay, information *info)
{
    std::unique_ptr<PQSetting_3D> actAllDomain = std::make_unique<PQSetting_3D>(info->DIMENSION);
    std::unique_ptr<PQSetting_3D> actGlobalInLocal = std::make_unique<PQSetting_3D>(info->DIMENSION);
    
    // physical quantity of global and local patch
    IIM_3D::makeDisplacementGradientAtGP(integralElementID, integralPointNaturalTilde, actAllDomain.get(), info);  // ガウス点での変位勾配を計算
    IIM_3D::makeStrainAtGP(integralElementID, integralPointNaturalTilde, actAllDomain.get(), info);                // ガウス点でのひずみを計算
    IIM_3D::makeStress(actAllDomain.get(), info);                                                                  // ガウス点での応力を計算

    // calculate the global patch physical quantity in local coordinate
    IIM_3D::checkGlobalInLocalCoordinates(integralElementID, gaussID, actGlobalInLocal.get(), info);               // ロー���ルガ���ス点がどのグローバル要素に属するかを確認
    IIM_3D::makeStress(actGlobalInLocal.get(), info);                                                              // ローカルガウス点座標に該当するグローバルでの応力を計算

    //overlay
    actOverlay->setDispGrad1(math_->overlay(actAllDomain->getDispGrad1(), actGlobalInLocal->getDispGrad1(), info->DIMENSION));  // ξ方向変位勾配の重ね合わせ
    actOverlay->setDispGrad2(math_->overlay(actAllDomain->getDispGrad2(), actGlobalInLocal->getDispGrad2(), info->DIMENSION));  // η方向変位勾配の重ね合わせ
    actOverlay->setDispGrad3(math_->overlay(actAllDomain->getDispGrad3(), actGlobalInLocal->getDispGrad3(), info->DIMENSION));  // ζ方向変位勾配の重ね合わせ
    actOverlay->setStrain(math_->overlay(actAllDomain->getStrain(), actGlobalInLocal->getStrain(), D_MATRIX_SIZE));             // ひずみの重ね合わせ
    actOverlay->setStress(math_->overlay(actAllDomain->getStress(), actGlobalInLocal->getStress(), D_MATRIX_SIZE));             // 応力の重ね合わせ
    IIM_3D::makeStrainEnergyDensityAtGP(actOverlay);                                                                            // 重ね場におけるひずみエネルギー密度の計算
}

void IIM_3D::computeAuxiliaryField(const int mode, vector<double> &JintegralCenter, vector<double> &integralPointCoords, AuxiliaryField *auxField, PQSetting_3D *actOverlay, PQSetting_3D *aux, information *info)
{
    int i;

    // orthogonal basis vector Q and Qt
    vector<double> Q(info->DIMENSION * info->DIMENSION, 0.0), Qt(info->DIMENSION * info->DIMENSION), tempTensor(info->DIMENSION * info->DIMENSION, 0.0);
    for (i = 0; i < info->DIMENSION; i++) 
    {
        Q[i * info->DIMENSION + 0] = crackFrontVector[2 * info->DIMENSION + i];  // 従法線基底ベクトル
        Q[i * info->DIMENSION + 1] = crackFrontVector[0 * info->DIMENSION + i];  // 法線基底ベクトル
        Q[i * info->DIMENSION + 2] = crackFrontVector[1 * info->DIMENSION + i];  // 接線基底ベクトル
        Qt[i * info->DIMENSION + i] = 1.0;                                       // 単位行列 [1 0 0; 0 1 0; 0 0 1]
    }

    // SIF計算点をx-y-z座標系からき裂方向座標 x'-y'-z'へと変換
    vector<double> localCenter(info->DIMENSION), localVertex(info->DIMENSION), localPoint(info->DIMENSION);  // local : x'-y'-z'座標系 [従法線-接線-法線]
    localCenter[0] = JintegralCenter[0] * crackFrontVector[2 * info->DIMENSION + 0] + JintegralCenter[1] * crackFrontVector[2 * info->DIMENSION + 1] + JintegralCenter[2] * crackFrontVector[2 * info->DIMENSION + 2];  // 従法線
    localCenter[1] = JintegralCenter[0] * crackFrontVector[0 * info->DIMENSION + 0] + JintegralCenter[1] * crackFrontVector[0 * info->DIMENSION + 1] + JintegralCenter[2] * crackFrontVector[0 * info->DIMENSION + 2];  // 法線
    localCenter[2] = JintegralCenter[0] * crackFrontVector[1 * info->DIMENSION + 0] + JintegralCenter[1] * crackFrontVector[1 * info->DIMENSION + 1] + JintegralCenter[2] * crackFrontVector[1 * info->DIMENSION + 2];  // 接線

    // 積分点をき裂方向座標へ変換
    localVertex[0] = integralPointCoords[0] * crackFrontVector[2 * info->DIMENSION + 0] + integralPointCoords[1] * crackFrontVector[2 * info->DIMENSION + 1] + integralPointCoords[2] * crackFrontVector[2 * info->DIMENSION + 2];
    localVertex[1] = integralPointCoords[0] * crackFrontVector[0 * info->DIMENSION + 0] + integralPointCoords[1] * crackFrontVector[0 * info->DIMENSION + 1] + integralPointCoords[2] * crackFrontVector[0 * info->DIMENSION + 2];
    localVertex[2] = integralPointCoords[0] * crackFrontVector[1 * info->DIMENSION + 0] + integralPointCoords[1] * crackFrontVector[1 * info->DIMENSION + 1] + integralPointCoords[2] * crackFrontVector[1 * info->DIMENSION + 2];

    for (i = 0; i < info->DIMENSION; i++) localPoint[i] = localVertex[i] - localCenter[i];  // 積分点からSIF計算点までのベクトル

    double R, theta;  
    vector<double> dR(info->DIMENSION), dTheta(info->DIMENSION);
    R = sqrt(localPoint[0] * localPoint[0] + localPoint[1] * localPoint[1]);  // 積分点からSIF計算点までの距離
    theta = atan2(localPoint[1], localPoint[0]);                              // 積分点からSIF計算点までの角度

    dR[0] = localPoint[0] / R;  // ∂R/∂x
    dR[1] = localPoint[1] / R;  // ∂R/∂y
    dR[2] = 0.0;                // ∂R/∂z

    dTheta[0] = -localPoint[1] / (R * R);  // ∂θ/∂x
    dTheta[1] =  localPoint[0] / (R * R);  // ∂θ/∂y
    dTheta[2] = 0.0;                       // ∂θ/∂z 

    // displacement gradient in x-y-z coordinates
    // makeDisplacementGradientAux1でx'-y'-z'座標における変位勾配を計算後，tensorTransformation関数を用いてx-y-z座標系に変換
    vector<double> tempDispGradient(info->DIMENSION * info->DIMENSION);
    switch (mode) {
        case 1:     // mode I
            tempDispGradient = math_->MathFormula::tensorTransformation(Q, Qt, auxField->AuxiliaryField::makeDisplacementGradientAux1(R, theta, dR, dTheta, info), info); 
            break;
        case 2:     // mode II
            tempDispGradient = math_->MathFormula::tensorTransformation(Q, Qt, auxField->AuxiliaryField::makeDisplacementGradientAux2(R, theta, dR, dTheta, info), info);
            break;
        case 3:     // mode III
            tempDispGradient = math_->MathFormula::tensorTransformation(Q, Qt, auxField->AuxiliaryField::makeDisplacementGradientAux3(R, theta, dR, dTheta, info), info);
            break;
    }
    // auxオブジェクトに仮想場の変位勾配を設定
    auxField->AuxiliaryField::setDisplacementGradientAux(tempDispGradient, aux, info);

    // strain in global coordinates
    auxField->AuxiliaryField::makeStrainAux(tempDispGradient, aux, info);

    // stress in global coordinates
    IIM_3D::makeStress(aux, info);

    // mutual strain energy density
    double W = 0.0;
    for (i = 0; i < N_STRAIN; i++) W += actOverlay->getStress()[i] * aux->getStrain()[i];  // W_mutual = σij(actual) * εij(auxiliary)
    aux->setStrainEnergyDensity(W);                                                        
}

// [0, 1]空間で与えられた座標点が特定の要素内に存在するかどうかを判定する関数
bool IIM_3D::isCoordinateInElement(const vector<double> &center, const int patchID, const int elementID, information *info)  
{   
    // 指定のノットスパン座標がなければ0を返す
    double residual = 1.0E-12;
    double knot_span = 0.0;

    // center[ξ, η, ζ] が該当要素の始点のノットスパンより大きくないか & 
    // center[ξ, η, ζ] が該当要素の終点のノットスパンより小さくないか &
    // ノットスパンが0でないか
    for (int i = 0; i < info->DIMENSION; i++) 
    {   
        knot_span = info->Position_Knots[info->Total_Knot_to_patch_dim[patchID * info->DIMENSION + i] + info->Order[patchID * info->DIMENSION + i] + info->ENC[elementID * info->DIMENSION + i] + 0] 
                         - info->Position_Knots[info->Total_Knot_to_patch_dim[patchID * info->DIMENSION + i] + info->Order[patchID * info->DIMENSION + i] + info->ENC[elementID * info->DIMENSION + i] + 1];
        
        if (!(center[i] >= (info->Position_Knots[info->Total_Knot_to_patch_dim[patchID * info->DIMENSION + i] + info->Order[patchID * info->DIMENSION + i] + info->ENC[elementID * info->DIMENSION + i] + 0] - residual) &&  
              center[i] < info->Position_Knots[info->Total_Knot_to_patch_dim[patchID * info->DIMENSION + i] + info->Order[patchID * info->DIMENSION + i] + info->ENC[elementID * info->DIMENSION + i] + 1] + residual ) && 
              fabs(knot_span) > residual) return false;    
    }
    return true;
}

void IIM_3D::makeBgradientMatrix(int elementID, vector<double> &BGradient1, vector<double> &BGradient2, vector<double> &BGradient3, vector<double> &localCoords, information *info)
{
    int i, j, k;
    vector<double> a(info->DIMENSION * info->DIMENSION, 0.0), b(info->DIMENSION * MAX_NO_CP_ON_ELEMENT, 0.0);

    for (i = 0; i < info->DIMENSION; i++)
        for (j = 0; j < info->DIMENSION; j++)
            for (k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[elementID]]; k++)
                a[i * info->DIMENSION + j] += dShape_func(k, j, localCoords.data(), elementID, info) * info->Node_Coordinate[info->Controlpoint_of_Element[elementID * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i];

    InverseMatrix_3x3(a.data());

    for (i = 0; i < info->DIMENSION; i++)
        for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[elementID]]; j++)
            for (k = 0; k < info->DIMENSION; k++)
                b[i * MAX_NO_CP_ON_ELEMENT + j] += a[k * info->DIMENSION + i] * dShape_func(j, k, localCoords.data(), elementID, info);

    #pragma omp parallel for if(OpenMP)
    for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[elementID]]; i++) 
    {   
        BGradient1[0 * MAX_KIEL_SIZE + (info->DIMENSION * i + 0)] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
        BGradient1[0 * MAX_KIEL_SIZE + (info->DIMENSION * i + 1)] = 0.0;
        BGradient1[0 * MAX_KIEL_SIZE + (info->DIMENSION * i + 2)] = 0.0;
        BGradient1[1 * MAX_KIEL_SIZE + (info->DIMENSION * i + 0)] = 0.0;
        BGradient1[1 * MAX_KIEL_SIZE + (info->DIMENSION * i + 1)] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
        BGradient1[1 * MAX_KIEL_SIZE + (info->DIMENSION * i + 2)] = 0.0;
        BGradient1[2 * MAX_KIEL_SIZE + (info->DIMENSION * i + 0)] = 0.0;
        BGradient1[2 * MAX_KIEL_SIZE + (info->DIMENSION * i + 1)] = 0.0;
        BGradient1[2 * MAX_KIEL_SIZE + (info->DIMENSION * i + 2)] = b[0 * MAX_NO_CP_ON_ELEMENT + i];

        BGradient2[0 * MAX_KIEL_SIZE + (info->DIMENSION * i + 0)] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
        BGradient2[0 * MAX_KIEL_SIZE + (info->DIMENSION * i + 1)] = 0.0;
        BGradient2[0 * MAX_KIEL_SIZE + (info->DIMENSION * i + 2)] = 0.0;
        BGradient2[1 * MAX_KIEL_SIZE + (info->DIMENSION * i + 0)] = 0.0;
        BGradient2[1 * MAX_KIEL_SIZE + (info->DIMENSION * i + 1)] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
        BGradient2[1 * MAX_KIEL_SIZE + (info->DIMENSION * i + 2)] = 0.0;
        BGradient2[2 * MAX_KIEL_SIZE + (info->DIMENSION * i + 0)] = 0.0;
        BGradient2[2 * MAX_KIEL_SIZE + (info->DIMENSION * i + 1)] = 0.0;
        BGradient2[2 * MAX_KIEL_SIZE + (info->DIMENSION * i + 2)] = b[1 * MAX_NO_CP_ON_ELEMENT + i];

        BGradient3[0 * MAX_KIEL_SIZE + (info->DIMENSION * i + 0)] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
        BGradient3[0 * MAX_KIEL_SIZE + (info->DIMENSION * i + 1)] = 0.0;
        BGradient3[0 * MAX_KIEL_SIZE + (info->DIMENSION * i + 2)] = 0.0;
        BGradient3[1 * MAX_KIEL_SIZE + (info->DIMENSION * i + 0)] = 0.0;
        BGradient3[1 * MAX_KIEL_SIZE + (info->DIMENSION * i + 1)] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
        BGradient3[1 * MAX_KIEL_SIZE + (info->DIMENSION * i + 2)] = 0.0;
        BGradient3[2 * MAX_KIEL_SIZE + (info->DIMENSION * i + 0)] = 0.0;
        BGradient3[2 * MAX_KIEL_SIZE + (info->DIMENSION * i + 1)] = 0.0;
        BGradient3[2 * MAX_KIEL_SIZE + (info->DIMENSION * i + 2)] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
    }
}

void IIM_3D::makeOrthogonalBasisVector(const vector<double> &vector1, const vector<double> &vector2, information *info)
{
    int i;
    double magnitude;  // ベクトルの大きさ(L2ノルム)
    vector<double> t(info->DIMENSION), n(info->DIMENSION), b(info->DIMENSION);

    // tangent vector(き裂前縁に対する接線基底ベクトル)
    for (i = 0; i < info->DIMENSION; i++) t[i] = vector1[i] + vector2[i];  // vector1と2の内積
    magnitude = sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
    for (i = 0; i < info->DIMENSION; i++) 
    {
        crackFrontVector[1 * info->DIMENSION + i] = t[i] / magnitude;
        // cout << "Crack Front Vector[1]" << "[" << i << "]: " << crackFrontVector[1 * info->DIMENSION + i] << endl;
    }

    // normal vector(き裂面に対する法線基底ベクトル)
    n = math_->crossProduct3x3(vector1, vector2, info);  // vector1と2の外積
    magnitude = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    for (i = 0; i < info->DIMENSION; i++) 
    {
        crackFrontVector[0 * info->DIMENSION + i] = n[i] / magnitude;
        // cout << "Crack Front Vector[0]" << "[" << i << "]: " << crackFrontVector[0 * info->DIMENSION + i] << endl;
    }

    // binormal vector(き裂前縁に対する法線基底ベクトル)
    b = math_->crossProduct3x3(n, t, info);
    magnitude = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
    for (i = 0; i < info->DIMENSION; i++) 
    {
        crackFrontVector[2 * info->DIMENSION + i] = b[i] / magnitude;
        // cout << "Crack Front Vector[2]" << "[" << i << "]: " << crackFrontVector[2 * info->DIMENSION + i] << endl;
    }
}

void IIM_3D::makeDisplacementGradientAtGP(const int elementID, vector<double> &localCoords, PQSetting_3D *act, information *info)
{
    int i, j;
    vector<double> U(MAX_KIEL_SIZE);
    vector<double> BGradient1(info->DIMENSION * MAX_KIEL_SIZE, 0.0), BGradient2(info->DIMENSION * MAX_KIEL_SIZE, 0.0), BGradient3(info->DIMENSION * MAX_KIEL_SIZE, 0.0);
    vector<double> tempDispGrad1(info->DIMENSION, 0.0), tempDispGrad2(info->DIMENSION, 0.0), tempDispGrad3(info->DIMENSION, 0.0);

    int  KIEL_SIZE = info->DIMENSION * info->No_Control_point_ON_ELEMENT[info->Element_patch[elementID]];

    #pragma omp parallel for if(OpenMP)
    for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[elementID]]; i++) {
        for (j = 0; j < info->DIMENSION; j++) {
            U[i * info->DIMENSION + j] = info->Displacement[info->Controlpoint_of_Element[elementID * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j];
        }
    }

    IIM_3D::makeBgradientMatrix(elementID, BGradient1, BGradient2, BGradient3, localCoords, info);
    // #pragma omp parallel for reduction(+:tempDispGrad1[:info->DIMENSION],tempDispGrad2[:info->DIMENSION],tempDispGrad3[:info->DIMENSION]) if(OpenMP)
    for (i = 0; i < info->DIMENSION; i++) {
        for (j = 0; j < KIEL_SIZE; j++) {
            tempDispGrad1[i] += BGradient1[i * MAX_KIEL_SIZE + j] * U[j];
            tempDispGrad2[i] += BGradient2[i * MAX_KIEL_SIZE + j] * U[j];
            tempDispGrad3[i] += BGradient3[i * MAX_KIEL_SIZE + j] * U[j];
        }
    }

    act->setDispGrad1(tempDispGrad1);
    act->setDispGrad2(tempDispGrad2);
    act->setDispGrad3(tempDispGrad3);
}

void IIM_3D::makeStrainAtGP(const int elementID, vector<double> &localCoords, PQSetting_3D *act, information *info)
{
    int i, j;
    vector<double> U(MAX_KIEL_SIZE);
    vector<double> B(D_MATRIX_SIZE * MAX_KIEL_SIZE);
    vector<double> tempStrain(D_MATRIX_SIZE, 0.0);
    
    #pragma omp parallel for if(OpenMP)
    for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[elementID]]; i++) 
    {
        for (j = 0; j < info->DIMENSION; j++) {
            U[i * info->DIMENSION + j] = info->Displacement[info->Controlpoint_of_Element[elementID * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j];
        }
    }


    int KIEL_SIZE = info->DIMENSION * info->No_Control_point_ON_ELEMENT[info->Element_patch[elementID]];
    Make_B_Matrix_anypoint(elementID, B.data(), localCoords.data(), info);
    // #pragma omp parallel for reduction(+:tempStrain[:D_MATRIX_SIZE]) if(OpenMP)
    for (i = 0; i < D_MATRIX_SIZE; i++) 
    {
        for (j = 0; j < KIEL_SIZE; j++) 
        {
            tempStrain[i] += B[i * MAX_KIEL_SIZE + j] * U[j];
        }
    }
    act->setStrain(tempStrain);
}

void IIM_3D::makeStress(PQSetting_3D *act, information *info)
{
    vector<double> tempStress(D_MATRIX_SIZE, 0.0);
    #pragma omp parallel for if(OpenMP)
    for (int i = 0; i < D_MATRIX_SIZE; i++) 
    {
        for (int j = 0; j < D_MATRIX_SIZE; j++) 
        {
            tempStress[i] += info->D[i * D_MATRIX_SIZE + j] * act->getStrain()[j];
        }
    }
    act->setStress(tempStress);
}

void IIM_3D::makeStrainEnergyDensityAtGP(PQSetting_3D *act)
{
    double tempStrainEnergyDensity = 0.0;
    #pragma omp parallel for reduction(+:tempStrainEnergyDensity) if(OpenMP)
    for (int i = 0; i < D_MATRIX_SIZE; i++) tempStrainEnergyDensity += (act->getStress()[i] * act->getStrain()[i]) / 2.0;
    act->setStrainEnergyDensity(tempStrainEnergyDensity);
}


void IIM_3D::checkGlobalInLocalCoordinates(const int elementID, const int gaussID, PQSetting_3D *globalInLoc, information *info)
{
    if (info->Element_mesh[elementID] > 0)  // ローカルパッチのみ処理の対処
    {
        const auto& overlapping_elements = info->eoi[elementID];
        
        if (!overlapping_elements.empty())  // ローカルメッシュに重なるグローバルパッチがあれば処理
        {
            for (int overlapping_element : overlapping_elements)  // ローカル要素に重なるグローバル要素分だけ処理
            {
                // ローカル要素ガウス点座標におけるグローバル要素の変位勾配とひずみを計算
                IIM_3D::makeDispGradientAndStrain(elementID, overlapping_element, gaussID, globalInLoc, info);
            }
        }
    }
}

void IIM_3D::makeDispGradientAndStrain(const int localElementID, const int globalElementID, const int gaussID, PQSetting_3D *globalInLoc, information *info)
{
    int i, j;
    int BDBJFlag = 0;

    vector<vector<double>> G_Gxi(MAX_POW_NG_EXTEND, vector<double>(MAX_DIMENSION));

	// determine if Local Gaussian points are included in the Global element
	// calculate physical coordinates of Local element at Gaussian point
	double data_result_shape[MAX_DIMENSION] = {0.0};
	for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[localElementID]]; i++) {
		double R_shape_func = Shape_func(i, &info->Gxi[gaussID * info->DIMENSION], localElementID, info);
		for (j = 0; j < info->DIMENSION; j++) {
			data_result_shape[j] += R_shape_func * info->Node_Coordinate[info->Controlpoint_of_Element[localElementID * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + j];
		}
	}

	// calculate parameter space coordinates of Local element at Gaussian points on the Global patch
	double para[MAX_DIMENSION] = {0.0};
	int patch_n = 0;
    if (info->DIMENSION == 2) {
        for (i = 0; i < info->Total_Patch_on_mesh[0]; i++) 
        {                   
            // int iteration = Calc_xi_eta(data_result_shape[0], data_result_shape[1], 
            //                             &info->Position_Knots[info->Total_Knot_to_patch_dim[i * info->DIMENSION + 0]], &info->Position_Knots[info->Total_Knot_to_patch_dim[i * info->DIMENSION + 1]],
            //                             info->No_Control_point[i * info->DIMENSION + 0], info->No_Control_point[i * info->DIMENSION + 1], 
            //                             info->Order[i * info->DIMENSION + 0], info->Order[i * info->DIMENSION + 1],
            //                             &para[0], &para[1], i, info);
            int iteration = calc_patch_parameter_coord(data_result_shape, i, para, info);
            patch_n = i;

            if (iteration == ERROR) {
                printf("Not converged\n");
                exit(1);
            }
        }
    } else if (info->DIMENSION == 3) {
        for (i = 0; i < info->Total_Patch_on_mesh[0]; i++) {
            // int iteration = Calc_xi_eta_zeta(data_result_shape[0], data_result_shape[1], data_result_shape[2],
            //                                  &info->Position_Knots[info->Total_Knot_to_patch_dim[i * info->DIMENSION + 0]], &info->Position_Knots[info->Total_Knot_to_patch_dim[i * info->DIMENSION + 1]], &info->Position_Knots[info->Total_Knot_to_patch_dim[i * info->DIMENSION + 2]],
            //                                  info->No_Control_point[i * info->DIMENSION + 0], info->No_Control_point[i * info->DIMENSION + 1], info->No_Control_point[i * info->DIMENSION + 2],
            //                                  info->Order[i * info->DIMENSION + 0], info->Order[i * info->DIMENSION + 1], info->Order[i * info->DIMENSION + 2],
            //                                  &para[0], &para[1], &para[2], i, info);
            int iteration = calc_patch_parameter_coord(data_result_shape, i, para, info);
            patch_n = i;

            if (iteration == ERROR) {
                printf("Not converged\n");
                exit(1);
            }
        }
    }

	//要素内外判定
    if (info->DIMENSION == 2) {
        if (para[0] >= info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 0] + info->Order[patch_n * info->DIMENSION + 0] + info->ENC[globalElementID * info->DIMENSION + 0] + 0] &&
            para[0] <  info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 0] + info->Order[patch_n * info->DIMENSION + 0] + info->ENC[globalElementID * info->DIMENSION + 0] + 1] &&
            para[1] >= info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 1] + info->Order[patch_n * info->DIMENSION + 1] + info->ENC[globalElementID * info->DIMENSION + 1] + 0] &&
            para[1] <  info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 1] + info->Order[patch_n * info->DIMENSION + 1] + info->ENC[globalElementID * info->DIMENSION + 1] + 1])	//要素内であるとき
        {
            BDBJFlag = 1;

            //親要素座標の算出
            G_Gxi[gaussID][0] = -1.0 + 2.0 * (para[0] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 0] + info->Order[patch_n * info->DIMENSION + 0] + info->ENC[globalElementID * info->DIMENSION + 0] + 0]) / (info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 0] + info->Order[patch_n * info->DIMENSION + 0] + info->ENC[globalElementID * info->DIMENSION + 0] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 0] + info->Order[patch_n * info->DIMENSION + 0] + info->ENC[globalElementID * info->DIMENSION + 0] + 0]);
            G_Gxi[gaussID][1] = -1.0 + 2.0 * (para[1] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 1] + info->Order[patch_n * info->DIMENSION + 1] + info->ENC[globalElementID * info->DIMENSION + 1] + 0]) / (info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 1] + info->Order[patch_n * info->DIMENSION + 1] + info->ENC[globalElementID * info->DIMENSION + 1] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 1] + info->Order[patch_n * info->DIMENSION + 1] + info->ENC[globalElementID * info->DIMENSION + 1] + 0]);
        }
        else	// 要素外であるとき
        {
            BDBJFlag = 0;
        }
    } else if (info->DIMENSION == 3) {
        if (para[0] >= info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 0] + info->Order[patch_n * info->DIMENSION + 0] + info->ENC[globalElementID * info->DIMENSION + 0]] &&
            para[0] <  info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 0] + info->Order[patch_n * info->DIMENSION + 0] + info->ENC[globalElementID * info->DIMENSION + 0] + 1] &&
            para[1] >= info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 1] + info->Order[patch_n * info->DIMENSION + 1] + info->ENC[globalElementID * info->DIMENSION + 1]] &&
            para[1] <  info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 1] + info->Order[patch_n * info->DIMENSION + 1] + info->ENC[globalElementID * info->DIMENSION + 1] + 1] &&
            para[2] >= info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 2] + info->Order[patch_n * info->DIMENSION + 2] + info->ENC[globalElementID * info->DIMENSION + 2]] &&
            para[2] <  info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 2] + info->Order[patch_n * info->DIMENSION + 2] + info->ENC[globalElementID * info->DIMENSION + 2] + 1])
        {
            BDBJFlag = 1;

            // 親要素座標の算出
            G_Gxi[gaussID][0] = - 1.0 + 2.0 * (para[0] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 0] + info->Order[patch_n * info->DIMENSION + 0] + info->ENC[globalElementID * info->DIMENSION + 0]]) / (info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 0] + info->Order[patch_n * info->DIMENSION + 0] + info->ENC[globalElementID * info->DIMENSION + 0] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 0] + info->Order[patch_n * info->DIMENSION + 0] + info->ENC[globalElementID * info->DIMENSION + 0]]);
            G_Gxi[gaussID][1] = - 1.0 + 2.0 * (para[1] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 1] + info->Order[patch_n * info->DIMENSION + 1] + info->ENC[globalElementID * info->DIMENSION + 1]]) / (info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 1] + info->Order[patch_n * info->DIMENSION + 1] + info->ENC[globalElementID * info->DIMENSION + 1] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 1] + info->Order[patch_n * info->DIMENSION + 1] + info->ENC[globalElementID * info->DIMENSION + 1]]);
            G_Gxi[gaussID][2] = - 1.0 + 2.0 * (para[2] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 2] + info->Order[patch_n * info->DIMENSION + 2] + info->ENC[globalElementID * info->DIMENSION + 2]]) / (info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 2] + info->Order[patch_n * info->DIMENSION + 2] + info->ENC[globalElementID * info->DIMENSION + 2] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + 2] + info->Order[patch_n * info->DIMENSION + 2] + info->ENC[globalElementID * info->DIMENSION + 2]]);
        }
        else
        {
            BDBJFlag = 0;
        }
    }


	// 要素内であるとき、次を計算
	if (BDBJFlag) {
        IIM_3D::makeDisplacementGradientAtGP(globalElementID, G_Gxi[gaussID], globalInLoc, info);
        IIM_3D::makeStrainAtGP(globalElementID, G_Gxi[gaussID], globalInLoc, info);
	}
}


// SIF情報のCSVファイルへの書き込み
void IIM_3D::writeSIFsToCSV(const string &fileName, const int Jcp, const vector<double> &calculatePointCoords, const double degree, information *info, CrackDataSetting_3D *crackData)
{
    static bool isFirstWrite = true;
    std::ofstream csvFile;
    
    // 関数呼び出し時にファイルの初期化とヘッダーの作成を行う
    if (isFirstWrite) {
        csvFile.open(fileName, std::ios::out | std::ios::trunc);
        if (!csvFile) {
            std::cout << "Error: cannot open file(" << fileName << ")" << std::endl;
            return;
        }
        
        // Write header based on MODE and domain types
        csvFile << "Jcp,x,y,z,degree";
        std::vector<std::string> modeLabels = {"K1", "K2", "K3"};
        for (int mode = 0; mode < MODE; ++mode) 
        {
            for (int domain = 0; domain < crackData->getdomain_types(); ++domain) 
            {
                csvFile << "," << modeLabels[mode] << "_all";
                if (domain > 0) {
                    csvFile << " - " << domain;
                }
            }
        }
        csvFile << std::endl;
        
        isFirstWrite = false;
    } else {
        csvFile.open(fileName, std::ios::out | std::ios::app);
        if (!csvFile) {
            std::cout << "Error: cannot open file(" << fileName << ")" << std::endl;
            return;
        }
    }

    // Write data
    csvFile << Jcp;
    for (int dim = 0; dim < info->DIMENSION; dim++) 
        csvFile << "," << calculatePointCoords[dim];
    csvFile << "," << degree;

    // Jcpに対応するMintegralSIFの値を書き込む
    for (int mode = 0; mode < MODE; ++mode) 
    {
        for (int domain = 0; domain < crackData->getdomain_types(); ++domain) 
        {
            int index = Jcp * (MODE * crackData->getdomain_types()) + mode * crackData->getdomain_types() + domain;
            csvFile << std::scientific << std::setprecision(15);
            csvFile << "," << getMintegralSIF()[index];
        }
    }
    csvFile << std::endl;

    csvFile.close();
}


void IIM_3D::readSIFsFromCSV(const std::string &fileName, CrackDataSetting_3D *crackData) 
{   
    size_t number_of_crack_data = 5;  // き裂前縁角度を含むデータ数

    // CSVファイルの読み込み
    std::ifstream file(fileName);
    if (!file) 
    {
        std::cerr << "Error: Cannot open file " << fileName << std::endl;
        exit(1);
    }

    std::string header;
    std::getline(file, header);  // ヘッダー行をスキップ

    std::string line;
    int point_count = 0;  // SIF計算点数
    
    // CSVファイルの各行を読みとり
    while (std::getline(file, line)) 
    {
        std::stringstream ss(line);
        std::string value;
        std::vector<double> row_data;

        // CSVの各列を解析
        while (std::getline(ss, value, ',')) 
        {
            row_data.push_back(std::stod(value));
        }

        // データ数の確認
        if (row_data.size() < number_of_crack_data + MODE * crackData->getdomain_types()) 
        {
            std::cerr << "Error: Invalid data format in line " << point_count + 1 << std::endl;
            exit(1);
        }

        // データの格納
        printf("\n=== SIF Calculation Point %d ===\n", point_count);
        printf("Point Number: %d\n", static_cast<int>(row_data[0]));
        printf("Point Number: %d\n", point_count);

        
        // 計算点座標の格納
        printf("Coordinates (x, y, z): %.6f, %.6f, %.6f\n", row_data[1], row_data[2], row_data[3]);
        // for (int dim = 0; dim < info->DIMENSION; dim++) 
        // {
        //     setCalculationPointCoords(row_data[dim + 1]);
        // }
        
        // き裂前縁角度の格納
        // setcrackFrontDegrees(row_data[4]);
        printf("Crack Front Angle: %.6f degrees\n", row_data[4]);

        // SIF値の格納
        for (int mode = 0; mode < MODE; mode++) 
        {
            printf("\nMode %d SIFs:\n", mode + 1);
            for (int domain = 0; domain < crackData->getdomain_types(); domain++) 
            {
                int index = point_count * (MODE * crackData->getdomain_types()) + mode * crackData->getdomain_types() + domain;
                setMintegralSIF(index, row_data[number_of_crack_data + mode * crackData->getdomain_types() + domain]);
                printf("Domain Size %d: %.6e\n", domain + 1, getMintegralSIF()[index]);
            }
        }
        point_count++;
    }
    std::cout << std::endl;
    std::cout << "Successfully read " << point_count << " SIF calculation points" << std::endl << std::endl;

    // 読み取ったSIFの数と計算点数が一致するか確認
    if (point_count != crackData->getJintegralPoints()) 
    {
        std::cerr << "Error: Number of SIF calculation points does not match" << std::endl;
        exit(1);
    }
}


vector<double> AuxiliaryField::makeDisplacementGradientAux1(const double r, const double theta, const vector<double> &dR, const vector<double> &dTheta, information *info)
{
    vector<double> dispGradient(info->DIMENSION * info->DIMENSION, 0.0);
    for (int i = 0; i < info->DIMENSION; i++) 
    {
        dispGradient[i * info->DIMENSION + 0] = dR[i] * ((1.0 + nu) / E) * (1.0 / sqrt(2.0 * M_PI * r)) * cos(theta / 2.0) * (1.0 - 2.0 * nu + sin(theta / 2.0) * sin(theta / 2.0)) + 
                                                dTheta[i] * ((1.0 + nu) / E) * sqrt(r / (2.0 * M_PI)) * (-sin(theta / 2.0)) * ((1.0 - 2.0 * nu + sin(theta / 2.0) * sin(theta / 2.0)) - 2.0 * cos(theta / 2.0) * cos(theta / 2.0));
        dispGradient[i * info->DIMENSION + 1] = dR[i] * ((1.0 + nu) / E) * (1.0 / sqrt(2.0 * M_PI * r)) * sin(theta / 2.0) * (2.0 - 2.0 * nu - cos(theta / 2.0) * cos(theta / 2.0)) + 
                                                dTheta[i] * ((1.0 + nu) / E) * sqrt(r / (2.0 * M_PI)) * cos(theta / 2.0) * ((2.0 - 2.0 * nu - cos(theta / 2.0) * cos(theta / 2.0)) + 2.0 * sin(theta / 2.0) * sin(theta / 2.0));
        dispGradient[i * info->DIMENSION + 2] = 0.0;
    }

    return dispGradient;
}

vector<double> AuxiliaryField::makeDisplacementGradientAux2(const double r, const double theta, const vector<double> &dR, const vector<double> &dTheta, information *info)
{
    vector<double> dispGradient(info->DIMENSION * info->DIMENSION, 0.0);
    for (int i = 0; i < info->DIMENSION; i++) 
    {
        dispGradient[i * info->DIMENSION + 0] = dR[i] * ((1.0 + nu) / E) * (1.0 / sqrt(2.0 * M_PI * r)) * sin(theta / 2.0) * (2.0 - 2.0 * nu + cos(theta / 2.0) * cos(theta / 2.0)) + 
                                                dTheta[i] * ((1.0 + nu) / E) * sqrt(r / (2.0 * M_PI)) * cos(theta / 2.0) * ((2.0 - 2.0 * nu + cos(theta / 2.0) * cos(theta / 2.0)) - 2.0 * sin(theta / 2.0) * sin(theta / 2.0));
        dispGradient[i * info->DIMENSION + 1] = dR[i] * ((1.0 + nu) / E) * (1.0 / sqrt(2.0 * M_PI * r)) * (-cos(theta / 2.0)) * (1.0 - 2.0 * nu - sin(theta / 2.0) * sin(theta / 2.0)) + 
                                                dTheta[i] * ((1.0 + nu) / E) * sqrt(r / (2.0 * M_PI)) * sin(theta / 2.0) * ((2.0 * nu + cos(theta / 2.0) * cos(theta / 2.0)) + 2.0 * cos(theta / 2.0) * cos(theta / 2.0));
        dispGradient[i * info->DIMENSION + 2] = 0.0;
    }

    return dispGradient;
}

vector<double> AuxiliaryField::makeDisplacementGradientAux3(const double r, const double theta, const vector<double> &dR, const vector<double> &dTheta, information *info)
{
    vector<double> dispGradient(info->DIMENSION * info->DIMENSION, 0.0);
    for (int i = 0; i < info->DIMENSION; i++) 
    {
        dispGradient[i * info->DIMENSION + 0] = 0.0;
        dispGradient[i * info->DIMENSION + 1] = 0.0;
        dispGradient[i * info->DIMENSION + 2] = dR[i] * 2.0 * ((1.0 + nu) / E) * (1.0 / sqrt(2.0 * M_PI * r)) * sin(theta / 2.0) +
                                                dTheta[i] * 2.0 * ((1.0 + nu) / E) * sqrt(r / (2.0 * M_PI)) * cos(theta / 2.0);
    }

    return dispGradient;
}

void AuxiliaryField::setDisplacementGradientAux(const vector<double> &tempDispGradient, PQSetting_3D *aux, information *info)
{
    vector<double> tempDispGradient1(info->DIMENSION), tempDispGradient2(info->DIMENSION), tempDispGradient3(info->DIMENSION);
    for (int i = 0; i < info->DIMENSION; i++) {
        tempDispGradient1[i] = tempDispGradient[0 * info->DIMENSION + i];
        tempDispGradient2[i] = tempDispGradient[1 * info->DIMENSION + i];
        tempDispGradient3[i] = tempDispGradient[2 * info->DIMENSION + i];
    }

    aux->setDispGrad1(tempDispGradient1);
    aux->setDispGrad2(tempDispGradient2);
    aux->setDispGrad3(tempDispGradient3);
}

void AuxiliaryField::makeStrainAux(const vector<double> &tempDispGradient, PQSetting_3D *aux, information *info)
{
    vector<double> tempStrain(info->DIMENSION * info->DIMENSION), strain(D_MATRIX_SIZE);

    #pragma omp parallel for if(OpenMP)
    for (int i = 0; i < info->DIMENSION; i++) 
    {
        for (int j = 0; j < info->DIMENSION; j++) 
        {
            if (i == j) {
                tempStrain[i * info->DIMENSION + j] = (tempDispGradient[i * info->DIMENSION + j] + tempDispGradient[j * info->DIMENSION + i]) / 2.0;
            } else {
                tempStrain[i * info->DIMENSION + j] = (tempDispGradient[i * info->DIMENSION + j] + tempDispGradient[j * info->DIMENSION + i]);
            }
        }
    }
    strain[0] = tempStrain[0 * info->DIMENSION + 0];    // xx
    strain[1] = tempStrain[1 * info->DIMENSION + 1];    // yy
    strain[2] = tempStrain[2 * info->DIMENSION + 2];    // zz
    strain[3] = tempStrain[0 * info->DIMENSION + 1];    // xy
    strain[4] = tempStrain[1 * info->DIMENSION + 2];    // yz
    strain[5] = tempStrain[2 * info->DIMENSION + 0];    // zx

    aux->setStrain(strain);
}


#if 1
/*
特異パッチの場合
q(ξ,η) = {
(1 - ξ / ξ1) * (1 - (η - η0) / (η1 - η0))  (η1 ≤ η < η0)
(1 - ξ / ξ1) * (1 - (η - η0) / (η2 - η0))  (η0 ≤ η < η2) }
をもとに∂q/∂x, ∂q/∂y, ∂q/∂zを計算

通常のパッチの場合
q(ξ,η,ζ) = {
(1 - ξ / ξ1) * (1 - ζ / ζ1) * (1 - (η - ζ0) * (1 - (η - η0) / (η1 - η0))   (η1 ≤ η < η0)
(1 - ξ / ξ1) * (1 - ζ / ζ1) * (1 - (η - ζ0) * (1 - (η - η0) / (η2 - η0))   (η0 ≤ η < η2) }
をもとに∂q/∂x, ∂q/∂y, ∂q/∂zを計算
*/
// q 関数の勾配を計算
void QGradientFunction::makeQGradient(int patchID, double crackJintegralNatural, const vector<double> &Jac, const vector<double> &dXiTilde_dXiPrime, const vector<double> &integralPointNatural, CrackDataSetting_3D *crackData, information *info)
{
    int i, domainType;
    vector<double> aJac;
    vector<double> dXi_dXiTilde(info->DIMENSION);                 // dξ' / dξ^~,  dη' / dη^~, dζ' / dζ^~
    vector<double> dq(info->DIMENSION);                           // dξ' / dξ^~, dη' / dη^~, dζ' / dζ^~ 
    vector<double> to_square_integralPointNatural(info->DIMENSION);  // き裂先端を原点とした座標系に変換した積分点座標

    std::copy(Jac.begin(), Jac.end(), std::back_inserter(aJac));  // aJac に Jac をコピー
    InverseMatrix_3x3(aJac.data());                               // aJac の逆行列

    for (i = 0; i < info->DIMENSION; i++)              
        dXi_dXiTilde[i] = 1.0 / dXiTilde_dXiPrime[i];  // dξ'/dξ^~ = 1 / (dξ^~ / dξ'), dη'/dη^~ = 1 / (dη^~ / dη'), dζ'/dζ^~ = 1 / (dζ^~ / dζ')
    
    // 領域サイズに応じた q 関数勾配を計算
    for (domainType = 0; domainType < crackData->getdomain_types(); domainType++)  
    {
        double domainXi = static_cast<double>(crackData->getIntegralDomainDivisions()[0] - domainType) / static_cast<double>(crackData->getIntegralDomainDivisions()[0]);    // ξ方向の領域サイズ[0 ~ 1]
        // double domainZeta = static_cast<double>(crackData->getIntegralDomainDivisions()[2] - domainType) / static_cast<double>(crackData->getIntegralDomainDivisions()[2]);  // ζ方向の領域サイズ[0 ~ 1] 
        double domainZeta = domainXi;  // ζ方向の領域サイズ[0 ~ 1] (ξ, η方向の要素数が同一であると仮定) 四辺形パッチの場合の計算にのみ使用  
        
        // ξ, ζ方向の領域サイズが負の値の場合エラー
        if (domainXi < 0.0 || domainXi > 1.0 || domainZeta < 0.0 || domainZeta > 1.0)
        {
            std::cout << "Error: domainXi or domainZeta is invalid value" << std::endl;
            std::cout << "domainXi = " << domainXi << ", domainZeta = " << domainZeta << std::endl;
            exit(1);
        }

        if (crackData->getpatchType() == 0)  
        {   
            if (integralPointNatural[1] >= crackJintegralNatural && integralPointNatural[0] <= domainXi)                                       // 積分点のη方向座標がSIF計算点以上の場合 (η_0 <= η < η_2), η_2 - η_0 = IntegralDomainNaturalWidth
            {                                     
                dq[0] = -(1.0 - ((integralPointNatural[1] - crackJintegralNatural) / crackData->getIntegralDomainNaturalWidth())) / domainXi;  // dq / dξ = -(1.0 - (η - η_0) / η_2 - η_0) / ξ_1    
                dq[1] = -(1.0 - integralPointNatural[0] / domainXi) / crackData->getIntegralDomainNaturalWidth();                              // dq / dη = -(1.0 - ξ / ξ_1) / (η_2 - η_0))
            } 
            else if (integralPointNatural[1] < crackJintegralNatural && integralPointNatural[0] <= domainXi)                                   // 積分点のη方向座標がSIF計算点より小さい場合 (η_1 <= η < η_0), η_1 - η_0 = - IntegralDomainNaturalWidth
            {                               
                dq[0] = -(1.0 + ((integralPointNatural[1] - crackJintegralNatural) / crackData->getIntegralDomainNaturalWidth())) / domainXi;  // dq / dξ = -(1.0 - (η - η_0) / η_1 - η_0) / ξ_1   
                dq[1] = (1.0 - integralPointNatural[0] / domainXi) / crackData->getIntegralDomainNaturalWidth();                               // dq / dη = -(1.0 - ξ / ξ_1) / (η_1 - η_0))
            } 
            else 
            {                                                                                                                                  // 設定積分領域範囲外はゼロ
                dq[0] = 0.0;                                                                                                                   
                dq[1] = 0.0;                                                                                                                   
            }
        } else  
        {   
            #if 1  // パッチ順パターン1
            if (patchID >= 1 && patchID <= crackData->getnum_edge_crack_patch())  // 第二象限のパッチ 
            {   
                to_square_integralPointNatural[0] = 1 - integralPointNatural[0];  // ξ ⇒ 1 - ξ
                to_square_integralPointNatural[1] = integralPointNatural[1];
                to_square_integralPointNatural[2] = integralPointNatural[2];
            }
            else if (patchID > crackData->getnum_edge_crack_patch() && patchID <= 2 * crackData->getnum_edge_crack_patch())  // 第一象限のパッチ
            {
                to_square_integralPointNatural[0] = integralPointNatural[0];  
                to_square_integralPointNatural[1] = integralPointNatural[1];
                to_square_integralPointNatural[2] = integralPointNatural[2];  
            }
            else if (patchID > 2 * crackData->getnum_edge_crack_patch() && patchID <= 3 * crackData->getnum_edge_crack_patch())  // 第三象限のパッチ
            {
                to_square_integralPointNatural[0] = 1 - integralPointNatural[0];  // ξ ⇒ 1 - ξ
                to_square_integralPointNatural[1] = integralPointNatural[1]; 
                to_square_integralPointNatural[2] = 1 - integralPointNatural[2];  // ζ ⇒ 1 - ζ
            }
            else if (patchID > 3 * crackData->getnum_edge_crack_patch() && patchID <= 4 * crackData->getnum_edge_crack_patch())  // 第四象限のパッチ
            {
                to_square_integralPointNatural[0] = integralPointNatural[0];  
                to_square_integralPointNatural[1] = integralPointNatural[1];
                to_square_integralPointNatural[2] = 1 - integralPointNatural[2];  // ζ ⇒ 1 - ζ
            }
            #endif 

            #if 0  // パッチ順パターン2
            if (patchID >= 1 && patchID <= crackData->getnum_edge_crack_patch())  // 第一象限のパッチ 
            {   
                to_square_integralPointNatural[0] = integralPointNatural[0];  
                to_square_integralPointNatural[1] = integralPointNatural[1];
                to_square_integralPointNatural[2] = integralPointNatural[2];
            }
            else if (patchID > crackData->getnum_edge_crack_patch() && patchID <= 2 * crackData->getnum_edge_crack_patch())  // 第二象限のパッチ
            {
                to_square_integralPointNatural[0] = 1 - integralPointNatural[0];  
                to_square_integralPointNatural[1] = integralPointNatural[1];
                to_square_integralPointNatural[2] = integralPointNatural[2];  
            }
            else if (patchID > 2 * crackData->getnum_edge_crack_patch() && patchID <= 3 * crackData->getnum_edge_crack_patch())  // 第四象限のパッチ
            {
                to_square_integralPointNatural[0] = integralPointNatural[0];  
                to_square_integralPointNatural[1] = integralPointNatural[1]; 
                to_square_integralPointNatural[2] = 1 - integralPointNatural[2];  // ζ ⇒ 1 - ζ
            }
            else if (patchID > 3 * crackData->getnum_edge_crack_patch() && patchID <= 4 * crackData->getnum_edge_crack_patch())  // 第三象限のパッチ
            {
                to_square_integralPointNatural[0] = 1 - integralPointNatural[0];  // ξ ⇒ 1 - ξ
                to_square_integralPointNatural[1] = integralPointNatural[1];
                to_square_integralPointNatural[2] = 1 - integralPointNatural[2];  // ζ ⇒ 1 - ζ
            }
            #endif 

            // dq/dξ, dq/dη, dq/dζ の計算
            if (to_square_integralPointNatural[1] >= crackJintegralNatural && to_square_integralPointNatural[0] <= domainXi && to_square_integralPointNatural[2] <= domainZeta)  // 積分点のη方向座標がSIF計算点以上の場合 (η_0 <= η < η_2), η_2 - η_0 = IntegralDomainNaturalWidth かつ ξ,ηの指定領域内の場合
            {           
                dq[0] = (-1 / domainXi) * 
                        (1 - to_square_integralPointNatural[2] / domainZeta) * 
                        (1.0 - (to_square_integralPointNatural[1] - crackJintegralNatural)
                        / crackData->getIntegralDomainNaturalWidth());                               // dq / dξ = - (1 - ζ / ζ_1) * (1.0 - (η - η_0) / η_2 - η_0) / ξ_1

                dq[1] = (1 - to_square_integralPointNatural[0] / domainXi) * 
                        (1 - to_square_integralPointNatural[2] / domainZeta) * 
                        (-1 / crackData->getIntegralDomainNaturalWidth());                     // dq / dη = -(1.0 - ξ / ξ_1) / (η_2 - η_0))
                
                dq[2] = (1 - to_square_integralPointNatural[0] / domainXi) * 
                        (-1  /domainZeta) * 
                        (1.0 - (to_square_integralPointNatural[1] - crackJintegralNatural)
                        / crackData->getIntegralDomainNaturalWidth());                       // dq / dζ = - (1 - ξ / ξ_1) * (1.0 - (η - η_0) / η_2 - η_0) / ζ_1
            } 
            else if (to_square_integralPointNatural[1] < crackJintegralNatural && to_square_integralPointNatural[0] <= domainXi && to_square_integralPointNatural[2] <= domainZeta)  // 積分点のη方向座標がSIF計算点より小さい場合 (η_1 <= η < η_0), η_1 - η_0 = - IntegralDomainNaturalWidth かつ ξ,ηの指定領域内の場合
            {                               
                dq[0] = (-1 / domainXi) * 
                        (1 - to_square_integralPointNatural[2] / domainZeta) * 
                        (1.0 + (to_square_integralPointNatural[1] - crackJintegralNatural)
                / crackData->getIntegralDomainNaturalWidth());                               // dq / dξ = - (1 - ζ / ζ_1) * (1.0 - (η - η_0) / η_2 - η_0) / ξ_1

                dq[1] = (1 - to_square_integralPointNatural[0] / domainXi) * 
                        (1 - to_square_integralPointNatural[2] / domainZeta) * 
                        (1 / crackData->getIntegralDomainNaturalWidth());                     // dq / dη = (1.0 - ξ / ξ_1) / (η_2 - η_0))
                
                dq[2] = (1 - to_square_integralPointNatural[0] / domainXi) * 
                        (-1 / domainZeta) * 
                        (1.0 + (to_square_integralPointNatural[1] - crackJintegralNatural)
                        / crackData->getIntegralDomainNaturalWidth());                       // dq / dζ = - (1 - ξ / ξ_1) * (1.0 - (η - η_0) / η_2 - η_0) / ζ_1                                         
            } 
            else   // 設定積分領域範囲外はゼロ
            {                                                                                                                           
                dq[0] = 0.0;                                                                                                                   
                dq[1] = 0.0;
                dq[2] = 0.0;                                                                                                                   
            }
        }


        for (i = 0; i < info->DIMENSION; i++) 
        {   
            if (crackData->getpatchType() == 0)
            {
                qGradient1[domainType * info->DIMENSION + i] = dq[0] * dXi_dXiTilde[0] * aJac[0 * info->DIMENSION + 0] + dq[1] * dXi_dXiTilde[1] * aJac[1 * info->DIMENSION + 0];  // ∂q/∂x = 
                qGradient2[domainType * info->DIMENSION + i] = dq[0] * dXi_dXiTilde[0] * aJac[0 * info->DIMENSION + 1] + dq[1] * dXi_dXiTilde[1] * aJac[1 * info->DIMENSION + 1];  // ∂q/∂y = 
                qGradient3[domainType * info->DIMENSION + i] = dq[0] * dXi_dXiTilde[0] * aJac[0 * info->DIMENSION + 2] + dq[1] * dXi_dXiTilde[1] * aJac[1 * info->DIMENSION + 2];  // ∂q/∂z = 
                // std::cout << "patchID : " << patchID << ", (ξ, η, ζ) : (" << integralPointNatural[0] << ", "<< integralPointNatural[1] << ", " << integralPointNatural[2] << ")" << ", qGradient[1] : " << qGradient1[domainType * info->DIMENSION + i] << std::endl;
                // std::cout << "patchID : " << patchID << ", (ξ, η, ζ) : (" << integralPointNatural[0] << ", "<< integralPointNatural[1] << ", " << integralPointNatural[2] << ")" << ", qGradient[2] : " << qGradient2[domainType * info->DIMENSION + i] << std::endl;
                // std::cout << "patchID : " << patchID << ", (ξ, η, ζ) : (" << integralPointNatural[0] << ", "<< integralPointNatural[1] << ", " << integralPointNatural[2] << ")" << ", qGradient[3] : " << qGradient3[domainType * info->DIMENSION + i] << std::endl;
            }
            else
            {   
                qGradient1[domainType * info->DIMENSION + i] = dq[0] * dXi_dXiTilde[0] * aJac[0 * info->DIMENSION + 0] + dq[1] * dXi_dXiTilde[1] * aJac[1 * info->DIMENSION + 0] + dq[2] * dXi_dXiTilde[2] * aJac[2 * info->DIMENSION + 0];  // ∂q/∂x = 
                qGradient2[domainType * info->DIMENSION + i] = dq[0] * dXi_dXiTilde[0] * aJac[0 * info->DIMENSION + 1] + dq[1] * dXi_dXiTilde[1] * aJac[1 * info->DIMENSION + 1] + dq[2] * dXi_dXiTilde[2] * aJac[2 * info->DIMENSION + 1];  // ∂q/∂y = 
                qGradient3[domainType * info->DIMENSION + i] = dq[0] * dXi_dXiTilde[0] * aJac[0 * info->DIMENSION + 2] + dq[1] * dXi_dXiTilde[1] * aJac[1 * info->DIMENSION + 2] + dq[2] * dXi_dXiTilde[2] * aJac[2 * info->DIMENSION + 2];  // ∂q/∂z = 
                // std::cout << "patchID : " << patchID << ", (ξ, η, ζ) : (" << to_square_integralPointNatural[0] << ", "<< to_square_integralPointNatural[1] << ", " << to_square_integralPointNatural[2] << ")" << ", qGradient[1] : " << qGradient1[domainType * info->DIMENSION + i] << std::endl;
                // std::cout << "patchID : " << patchID << ", (ξ, η, ζ) : (" << to_square_integralPointNatural[0] << ", "<< to_square_integralPointNatural[1] << ", " << to_square_integralPointNatural[2] << ")" << ", qGradient[2] : " << qGradient2[domainType * info->DIMENSION + i] << std::endl;
                // std::cout << "patchID : " << patchID << ", (ξ, η, ζ) : (" << to_square_integralPointNatural[0] << ", "<< to_square_integralPointNatural[1] << ", " << to_square_integralPointNatural[2] << ")" << ", qGradient[3] : " << qGradient3[domainType * info->DIMENSION + i] << std::endl;
            }
        }
    }
}
#endif 


// 同一要素内の[0, 1]空間で与えられた座標点0を[-1, 1]座標へ変換
vector<double> MathFormula::convertNaturalCoordinatesToTilde(const int patchID, const int elementID, const vector<double> &naturalCoords, information *info)
{
    vector<double> tildeCoords(info->DIMENSION, 0.0);
    for (int i = 0; i < info->DIMENSION; i++) 
    {
        tildeCoords[i] = 2.0 * (naturalCoords[i] - info->Position_Knots[info->Total_Knot_to_patch_dim[patchID * info->DIMENSION + i] + info->Order[patchID * info->DIMENSION + i] + info->ENC[elementID * info->DIMENSION + i] + 0]) / (info->Position_Knots[info->Total_Knot_to_patch_dim[patchID * info->DIMENSION + i] + info->Order[patchID * info->DIMENSION + i] + info->ENC[elementID * info->DIMENSION + i] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[patchID * info->DIMENSION + i] + info->Order[patchID * info->DIMENSION + i] + info->ENC[elementID * info->DIMENSION + i] + 0]) - 1.0;
    }
    // cout << "tildeCoords : " << tildeCoords[0] << ", " << tildeCoords[1] << ", " << tildeCoords[2] << endl;
    return tildeCoords;
}


// 3×3の外積
vector<double> MathFormula::crossProduct3x3(const vector<double> &a, const vector<double> &b, information *info)
{
    vector<double> c(info->DIMENSION);

    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];

    return c;
}


// き裂前縁x'-y'-z'座標系をx-y-z座標系へ変換
vector<double> MathFormula::tensorTransformation(const vector<double> &unitLocal, const vector<double> &unitGlobal, const vector<double> &inputTensor, information *info)
{
    int i, j, k;
    vector<double> Q(info->DIMENSION * info->DIMENSION, 0.0), Qt(info->DIMENSION * info->DIMENSION, 0.0), tempTensor(info->DIMENSION * info->DIMENSION, 0.0), outputTensor(info->DIMENSION * info->DIMENSION, 0.0);

    // 座標変換行列 Q[i][j] = Σ(unit_x'y'z'[k][i] * unit_xyz[k][j])
    for (i = 0; i < info->DIMENSION; i++)
        for (j = 0; j < info->DIMENSION; j++)
            for (k = 0; k < info->DIMENSION; k++)
                Q[j * info->DIMENSION + i] += unitLocal[k * info->DIMENSION + i] * unitGlobal[k * info->DIMENSION + j];

    // Qt[i][j] = Q[j][i] = Q^T
    for (i = 0; i < info->DIMENSION; i++)
        for (j = 0; j < info->DIMENSION; j++)
            Qt[i * info->DIMENSION + j] = Q[j * info->DIMENSION + i];
    
    // tempTensor_ij = Σ_k (Qt_ik * inputTensor_kj) = Q^T * inputTensor
    for (i = 0; i < info->DIMENSION; i++)
        for (j = 0; j < info->DIMENSION; j++)
            for (k = 0; k < info->DIMENSION; k++)
                tempTensor[j * info->DIMENSION + i] += Qt[k * info->DIMENSION + i] * inputTensor[j * info->DIMENSION + k];

    // outputTensor_ij = Σ_k (tempTensor_ik * Q_kj) = Q^T * inputTensor * Q
    for (i = 0; i < info->DIMENSION; i++)
        for (j = 0; j < info->DIMENSION; j++)
            for (k = 0; k < info->DIMENSION; k++)
                outputTensor[j * info->DIMENSION + i] += tempTensor[k * info->DIMENSION + i] * Q[j * info->DIMENSION + k];
    
    return outputTensor;
}


// 
vector<double> MathFormula::overlay(const vector<double> &value1, const vector<double> &value2, const int num)
{
    vector<double> value3(num);
    for (int i = 0; i < num; i++) value3[i] = value1[i] + value2[i];
    return value3;
}


