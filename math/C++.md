# 括号内嵌深度

## 问题描述

当一段字符存在小括号，中括号，大括号，判定括号的最大内嵌深度，如"(2+3)*[6+(1+1)/(1+1)]"最大深度为2。当表达式不合法时，如"]["返回0。

## 求解思路

利用“栈”的特点，合法的表达式必然是左括号"(/[/{"按顺序入栈，对应的右括号")/]/}"按顺序出栈。

```C++
#include<iostream>
#include<stack>
#include<algorithm>
using namespace std;

stack<char> myStack;
string str;
int flag=1;
int maxDepth=0;
void calDepth(int i){
    if(i==str.length()){
        return;
    }
    switch(str[i]){
        case '(':myStack.push(')');break;
        case '[':myStack.push(']');break;
        case '{':myStack.push('}');break;
        case ')':
        case ']':
        case '}':
                if(myStack.empty()||str[i]!=myStack.top()){
                    flag=0;
                    return;
                }
                else{
                    int depth=myStack.size();
                    maxDepth=max(maxDepth,depth);
                    myStack.pop();
                }
                break;
        }    
    calDepth(i+1);
}
int main(int args,char* argv[]){
    while(getline(cin,str)){
        calDepth(0);
        if(flag){
            cout<<maxDepth<<endl;
        }else{
            cout<<'0'<<endl; 
        }
        myStack.empty();
        maxDepth=0;
        flag=1;
    }
    return 0;
}
```

