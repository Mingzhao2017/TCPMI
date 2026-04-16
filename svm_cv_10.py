import sklearn
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, average_precision_score
from sklearn.svm import SVC
import pandas as pd
'''
data_names = {
    'DNA4mC_A_thaliana_3956', 'DNA4mC_C_elegans_3108', 'DNA4mC_C_equisetifolia_366',
    'DNA4mC_D_melanogaster_3538', 'DNA4mC_E_coli_1_776', 'DNA4mC_E_coli_2_3882',
    'DNA4mC_E_coli_3816', 'DNA4mC_F_vesca_15798', 'DNA4mC_G_pickeringii_1_1138',
    'DNA4mC_G_pickeringii_2_9028', 'DNA4mC_G_pickeringii_7522', 'DNA4mC_G_subterraneus_14128',
    'DNA4mC_G_subterraneus_1_1810', 'DNA4mC_G_subterraneus_2_19868', 'DNA4mC_S_cerevisiae_1980',
    'DNA4mC_Tolypocladium_15328', 'DNA6mA_C_elegans_7962', 'DNA6mA_C_equisetifolia_6066',
    'DNA6mA_D_melanogaster_11192', 'DNA6mA_F_vesca_1_3102', 'DNA6mA_F_vesca_8575',
    'DNA6mA_H_sapiens_18336', 'DNA6mA_R_chinensis_1_600', 'DNA6mA_R_chinensis_2134',
    'DNA6mA_S_cerevisiae_3786', 'DNA6mA_Tolypocladium_3380', 'DNA6mA_Xoc_BLS256_17216'
}
data_names = {'h_b_all', 'h_k_all', 'h_l_all', 'm_b_all', 'm_h_all', 'm_k_all', 'm_l_all', 'm_t_all', 'r_b_all', 'r_k_all', 'r_l_all'
}
'''
data_names = {'Lv'}
all_results = []
for dataname in data_names:
    train_matrix = np.loadtxt('./result/data_' + dataname + '_Tri_BiPSDP.csv', delimiter=' ')
    train_label = np.loadtxt('./result/label_' + dataname + '.csv', delimiter=' ')
    print(train_matrix.shape, train_label.shape)

    # 初始化结果列表
    results = []
    for i in range(0, 10):
        idx = 0
        cv_train_matrix = []
        cv_train_label = []
        cv_test_matrix = []
        cv_test_label = []
        idt = 0
        for line in train_matrix:
            if idx % 10 == i:
                cv_test_matrix.append(train_matrix[idt])
                cv_test_label.append(train_label[idt])
            else:
                cv_train_matrix.append(train_matrix[idt])
                cv_train_label.append(train_label[idt])
            idt += 1
            idx += 1

        # 使用SVM分类器
        svm_classifier = SVC(kernel='rbf', probability=True, C=0.25)  # 你可以根据需要调整参数
        svm_classifier.fit(cv_train_matrix, cv_train_label)
        cv_test_pred = svm_classifier.predict(cv_test_matrix)

        accuracy = accuracy_score(cv_test_label, cv_test_pred)
        print(i,accuracy)
        # 计算混淆矩阵
        tn, fp, fn, tp = confusion_matrix(cv_test_label, cv_test_pred).ravel()

        # 计算灵敏度（Sensitivity）
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0

        # 计算特异性（Specificity）
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0

        # 计算马修斯相关系数（MCC）
        mcc_value = (tp * tn - fp * fn) / np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) if (tp + fp) * (
                    tp + fn) * (tn + fp) * (tn + fn) > 0 else 0

        # 计算AUROC
        auroc_value = roc_auc_score(cv_test_label, svm_classifier.predict_proba(cv_test_matrix)[:, 1])

        # 计算AUPRC
        auprc_value = average_precision_score(cv_test_label, svm_classifier.predict_proba(cv_test_matrix)[:, 1])

        # 将结果添加到列表中
        results.append({
            'Fold': i,
            'Accuracy': accuracy,
            'Sensitivity': sensitivity,
            'Specificity': specificity,
            'MCC': mcc_value,
            'AUROC': auroc_value,
            'AUPRC': auprc_value
        })

    avg_results = {key: np.mean([result[key] for result in results]) for key in results[0] if key != 'Fold'}
    avg_results['Dataset'] = dataname  # 添加数据集名称
    all_results.append(avg_results)

# 创建 DataFrame
all_results_df = pd.DataFrame(all_results)

# 保存到 Excel 文件
all_results_df.to_excel('average_model_results_svm_0.5_Lv.xlsx', index=False)