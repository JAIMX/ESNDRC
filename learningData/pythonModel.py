
# coding: utf-8

# In[20]:


import pandas as pd


# In[21]:


data=pd.read_csv("training_small.csv")
data=data.drop_duplicates()


# In[56]:


data.head()


# In[26]:


# from sklearn.model_selection import train_test_split
# train_data,test_data=train_test_split(data,train_size=0.67,random_state=0)

# y_train=train_data.y
# train_data=train_data.drop(['y'],axis=1)
# x_train=train_data.as_matrix()

# y_test=test_data.y
# test_data=test_data.drop(['y'],axis=1)
# x_test=test_data.as_matrix()
# print(train_data.columns.tolist())


# In[88]:


from sklearn.linear_model import LogisticRegressionCV
y=data.y
x=data.as_matrix()
x=x[:,0:x.shape[1]-1]

clf=LogisticRegressionCV(cv=10,random_state=0,penalty="l1",solver="liblinear").fit(x,y)

coef_list=clf.coef_
col_names=data.columns.values.tolist()

for i in range(len(coef_list[0])):
    print(col_names[i])
    print(coef_list[0][i])
    print()

