#include <iostream>
#include "config.h"
#include "cpp_header.h"
#include <unordered_map>

using namespace std;

SMatF* make_binary(SMatF* mat)
{
	SMatF* one_mat = new SMatF(mat);	
	for(_int col=0;col<one_mat->nc;col++)
	{
		for(_int i=0;i<one_mat->size[col];i++)
			one_mat->data[col][i].second = 1;
	}
	return one_mat;
}

_float f(_float x,_float alpha)
{
	if(x==1)
		_float x = x-pow(10,-10);			
	return pow(x,alpha)*log(x);
}

SMatF* transform_matf(SMatF* mat,_float alpha)
{
	SMatF* ans = new SMatF(mat);
	for(_int col=0; col<mat->nc; col++)	
		for(_int i=0;i<mat->size[col];i++)	
			ans->data[col][i].second = f(mat->data[col][i].second, alpha);
	return ans;
}

_float g(_float x,_float alpha)
{
	return pow(x,alpha);
}

SMatF* transform_matg(SMatF* mat,_float alpha)
{
	SMatF* ans = new SMatF(mat);
	for(_int col=0; col<mat->nc; col++)	
		for(_int i=0;i<mat->size[col];i++)	
			ans->data[col][i].second = g(mat->data[col][i].second, alpha);
	return ans;
}

_float getvalue(int I, int J, SMatF* smat)
{
    for(int it = 0; it < smat->size[J]; ++it)
        if(smat->data[J][it].first == I)
            return smat->data[J][it].second;        
    return 0;
}

_int find_ind(int I, int J, SMatF* smat, _int start_ind)
{
    for(int it = start_ind; it < smat->size[J]; ++it)
    {
        if(smat->data[J][it].first == I)
            return it;        
        if(smat->data[J][it].first > I)
        	return -1;
    }
    return -1;
}

SMatF* transform_log(SMatF* mat)
{
	SMatF* ans = new SMatF(mat);
	for(_int col=0; col<mat->nc; col++)	
		for(_int i=0;i<mat->size[col];i++)
		{
			_float val = mat->data[col][i].second;
			if(abs(val-1)<1e-5)
				val = 1-1e-5;
			ans->data[col][i].second = log(val);
		}				
	return ans;
}

SMatF* transform_exp(SMatF* mat)
{
	SMatF* ans = new SMatF(mat);
	for(_int col=0; col<mat->nc; col++)	
		for(_int i=0;i<mat->size[col];i++)	
			ans->data[col][i].second = exp(mat->data[col][i].second);
	return ans;
}

SMatF* my_GM( SMatF* score_mat, SMatF* Y_token, _float alpha )    // GM new
{
	cout << alpha << "\nGM + alpha * logn"<<endl;
	
	SMatF* Y_token_t = Y_token->transpose();	
	SMatF* log_score_mat = transform_log(score_mat);	
	cout << "reached here"<<endl;
	_int dps = score_mat -> nc;
	_int labs = Y_token -> nc;

	SMatF* final_mat = new SMatF(labs, dps);

	for(_int dp=0;dp<dps;dp++)
	{				
		if(dp%1000==0)
			cout << dp << endl;
		unordered_map<_int,_int> umap;
		unordered_map<_int,_float> log_tok_scores;		

		unordered_map<_int,_int>::iterator itr;

		_int num_toks = score_mat->size[dp];
		for(_int i=0;i<num_toks;i++)
		{			
			_int tok = score_mat->data[dp][i].first;			
			log_tok_scores[tok] = log_score_mat->data[dp][i].second;			
			
			for(_int j=0; j< Y_token_t -> size[tok]; j++)
			{
				_int lab = Y_token_t->data[tok][j].first;
				if (umap.find(lab) == umap.end())
					umap[lab] = 1;
				else
					umap[lab]++;
			}
		}
		
		vector<int> curr_labs;
		for (itr = umap.begin(); itr != umap.end(); itr++) 
    	{
    		if(itr->second == Y_token -> size[itr->first])
    			curr_labs.push_back(itr->first);
    	} 
    	
		vector<pairIF> cur_col;
		for(_int i=0;i<curr_labs.size();i++)
		{		
			_int lab = curr_labs[i];

			_int toks = Y_token->size[lab];	

			pair<_int,_float>* vec_lab = Y_token->data[lab];	

			_float prodlog = 0;			

			for(_int j=0;j<toks;j++)
			{								
				_int tok = vec_lab[j].first;
				
				prodlog += vec_lab[j].second * log_tok_scores[tok];				
			}

			if(toks!=0)                                       // +clogn term
				prodlog += alpha*log(toks);

			cur_col.push_back(pairIF(lab, prodlog ));
		}

		sort(cur_col.begin(), cur_col.end(), comp_pair_by_second_desc<_int,_float>);
		cur_col.resize(min((size_t)100, cur_col.size()));

		final_mat->size[dp] = cur_col.size();
		final_mat->data[dp] = new pairIF[cur_col.size()];

		for(int i = 0; i < cur_col.size(); ++i)
			final_mat->data[dp][i] = cur_col[i];
	}

	SMatF* ans = transform_exp(final_mat);
	return ans;
}

//LABEL_LENGTH*100 WORKING

SMatF* my_func( SMatF* score_mat, SMatF* Y_token, _float alpha )
{
	cout << alpha << " x^-alpha.log(x)"<<endl;
	
	Y_token = make_binary(Y_token); 
	SMatF* Y_token_t = Y_token->transpose();	
	SMatF* f_score_mat = transform_matf(score_mat,alpha);
	SMatF* g_score_mat = transform_matg(score_mat,alpha);
	cout << "reached here"<<endl;
	_int dps = score_mat -> nc;
	_int labs = Y_token -> nc;

	SMatF* final_mat = new SMatF(labs, dps);

	for(_int dp=0;dp<dps;dp++)
	{				
		if(dp%1000==0)
			cout << dp << endl;
		unordered_map<_int,_int> umap;
		unordered_map<_int,_float> f_tok_scores;
		unordered_map<_int,_float> g_tok_scores;

		unordered_map<_int,_int>::iterator itr;

		_int num_toks = score_mat->size[dp];
		for(_int i=0;i<num_toks;i++)
		{			
			_int tok = score_mat->data[dp][i].first;			
			f_tok_scores[tok] = f_score_mat->data[dp][i].second;
			g_tok_scores[tok] = g_score_mat->data[dp][i].second;
			
			for(_int j=0; j< Y_token_t -> size[tok]; j++)
			{
				_int lab = Y_token_t->data[tok][j].first;
				if (umap.find(lab) == umap.end())
					umap[lab] = 1;
				else
					umap[lab]++;
			}
		}
		
		vector<int> curr_labs;
		for (itr = umap.begin(); itr != umap.end(); itr++) 
    	{
    		if(itr->second == Y_token -> size[itr->first])
    			curr_labs.push_back(itr->first);
    	} 
    	
		vector<pairIF> cur_col;
		for(_int i=0;i<curr_labs.size();i++)
		{		
			_int lab = curr_labs[i];

			_int toks = Y_token->size[lab];	

			pair<_int,_float>* vec_lab = Y_token->data[lab];	

			_float prodf = 0, prodg = 0;
			_int f_ind = 0, g_ind=0;

			for(_int j=0;j<toks;j++)
			{								
				_int tok = vec_lab[j].first;
				
				prodf += vec_lab[j].second * f_tok_scores[tok];
				prodg += vec_lab[j].second * g_tok_scores[tok];
			}
			//if(abs(prodf)<1e-10)
			//	cout << "prob=1\n";
			cur_col.push_back(pairIF(lab, prodf / prodg));
			
		}

		sort(cur_col.begin(), cur_col.end(), comp_pair_by_second_desc<_int,_float>);
		cur_col.resize(min((size_t)100, cur_col.size()));

		final_mat->size[dp] = cur_col.size();
		final_mat->data[dp] = new pairIF[cur_col.size()];

		for(int i = 0; i < cur_col.size(); ++i)
			final_mat->data[dp][i] = cur_col[i];
	}

	SMatF* ans = transform_exp(final_mat);
	return ans;
}

/*
this is L*100 DO NOT TOUCH
SMatF* my_func( SMatF* score_mat, SMatF* Y_token, _float alpha )
{
	cout << alpha <<endl;
	
	Y_token = make_binary(Y_token); 
	SMatF* Y_token_t = Y_token->transpose();	
	SMatF* f_score_mat = transform_matf(score_mat,alpha);
	SMatF* g_score_mat = transform_matg(score_mat,alpha);
	cout << "reached here"<<endl;
	_int dps = score_mat -> nc;
	_int labs = Y_token -> nc;

	SMatF* final_mat = new SMatF(labs, dps);

	for(_int dp=0;dp<dps;dp++)
	{		
		///////////////		
		if(dp%1000==0)
			cout << dp << endl;
		unordered_map<_int,_int> umap;
		unordered_map<_int,_int>::iterator itr;
		_int num_toks = score_mat->size[dp];
		for(_int i=0;i<num_toks;i++)
		{
			_int tok = score_mat->data[dp][i].first;
			//cout << "tok is "<<tok<<endl;
			for(_int j=0; j< Y_token_t -> size[tok]; j++)
			{
				_int lab = Y_token_t->data[tok][j].first;
				if (umap.find(lab) == umap.end())
					umap[lab] = 1;
				else
					umap[lab]++;
			}
		}
		/////////////
		vector<int> curr_labs;
		for (itr = umap.begin(); itr != umap.end(); itr++) 
    	{
    		if(itr->second == Y_token -> size[itr->first])
    			curr_labs.push_back(itr->first);
    	} 
    	//for(auto i = 0; i<curr_labs.size(); i++ )
    	//	cout << "currlabs "<<curr_labs[i]<<endl;
    	//cout << curr_labs << endl;
    	////////////
		vector<pairIF> cur_col;
		for(_int i=0;i<curr_labs.size();i++)
		{		
			_int lab = curr_labs[i];
		//for(_int lab=0;lab<labs;lab++)
		//{		
			_int toks = Y_token->size[lab];	
			pair<_int,_float>* vec_lab = Y_token->data[lab];	
			_float prodf = 0, prodg = 0;
			_int f_ind = 0, g_ind=0;
			for(_int tok=0;tok<toks;tok++)
			{
				f_ind = find_ind(vec_lab[tok].first,dp,f_score_mat,f_ind);				
				if(f_ind == -1)
				{
					cout << lab << "  " << tok << endl;				
					prodf=0;
					break;	
				}
				
				_float f_val = f_score_mat->data[dp][f_ind].second;
				prodf += vec_lab[tok].second * f_val;

				g_ind = find_ind(vec_lab[tok].first,dp,g_score_mat,g_ind);
				_float g_val = g_score_mat->data[dp][g_ind].second;
				prodg += vec_lab[tok].second * g_val;								
			}
			if(abs(prodf) > 1e-10)
				cur_col.push_back(pairIF(lab, prodf / prodg));
		}

		sort(cur_col.begin(), cur_col.end(), comp_pair_by_second_desc<_int,_float>);
		cur_col.resize(min((size_t)100, cur_col.size()));

		final_mat->size[dp] = cur_col.size();
		final_mat->data[dp] = new pairIF[cur_col.size()];

		for(int i = 0; i < cur_col.size(); ++i)
			final_mat->data[dp][i] = cur_col[i];
	}

	SMatF* ans = transform_exp(final_mat);
	return ans;
}




/*
SMatF* my_func( SMatF* score_mat, SMatF* Y_token, _float alpha )
{
	cout << alpha <<endl;
	
	Y_token = make_binary(Y_token); 
	SMatF* Y_token_t = Y_token->transpose();	
	SMatF* f_score_mat = transform_matf(score_mat,alpha);
	SMatF* g_score_mat = transform_matg(score_mat,alpha);
	cout << "reached here"<<endl;
	_int dps = score_mat -> nc;
	_int labs = Y_token -> nc;
	cout << dps << labs<<endl;
	SMatF* final_mat = new SMatF(labs, dps);
	
	vector<_int> v[labs];
	for(_int dp=0;dp<dps;dp++)
	{				
		_int num_toks = score_mat->size[dp];
		/*
		unordered_map<_int,_int> umap;
		unordered_map<_int,_int>::iterator itr;
		
		for(_int i=0;i<num_toks;i++)
		{
			_int tok = score_mat->data[dp][i].first;
			for(_int j=0; j< Y_token_t -> size[tok]; j++)
			{
				_int lab = Y_token_t->data[tok][j].first;
				if (umap.find(lab) == umap.end())
					umap[lab] = 1;
				else
					umap[lab]++;
			}
		}
		/*
		vector<int> curr_labs;
		for (itr = umap.begin(); itr != umap.end(); itr++) 
    	{
    		if(itr->second==num_toks)
    			curr_labs.push_back(itr->first);
    	} 
		*/
/*
		vector<pairIF> cur_col;
		for(_int lab=0;lab<labs;lab++)
		{			
		//for(_int i=0;i<curr_labs.size();i++)
		//{		
			//_int lab = curr_labs[i];
			_int toks = Y_token->size[lab];	
			pair<_int,_float>* vec_lab = Y_token->data[lab];	
			_float prodf = 0, prodg = 0;
			_int f_ind = 0, g_ind=0;
			for(_int tok=0;tok<toks;tok++)
			{
				f_ind = find_ind(vec_lab[tok].first,dp,f_score_mat,f_ind);	
				if(f_ind == -1)
				{
					prodf=0;
					break;	
				}
				_float f_val = f_score_mat->data[dp][f_ind].second;
				prodf += vec_lab[tok].second * f_val;

				g_ind = find_ind(vec_lab[tok].first,dp,g_score_mat,g_ind);
				_float g_val = g_score_mat->data[dp][g_ind].second;
				prodg += vec_lab[tok].second * g_val;								
			}
			if(abs(prodf) > 1e-10)
				cur_col.push_back(pairIF(lab, prodf / prodg));
		}

		sort(cur_col.begin(), cur_col.end(), comp_pair_by_second_desc<_int,_float>);
		cur_col.resize(min((size_t)100, cur_col.size()));

		final_mat->size[dp] = cur_col.size();
		final_mat->data[dp] = new pairIF[cur_col.size()];

		for(int i = 0; i < cur_col.size(); ++i)
			final_mat->data[dp][i] = cur_col[i];
		
	}
	cout << "here";
	SMatF* ans = transform_exp(final_mat);
	return ans;
}
*/
/*
SMatF* eliminate_zeros(SMatF* mat)
{
	_int cols = mat->nc;
	SMatF* final_mat = new SMatF(mat->nr,cols);
	for(_int col=0;col<cols;col++)
	{
		vector<pairIF> cur_col;
		for(_int i=0;i<mat->size[col];i++)
		{
			_int ind = mat->data[col][i].first;
			_float val = mat->data[col][i].second;
			if(abs(val)>1e-10)
				cur_col.push_back(pairIF(ind,val));
		}
		final_mat->size[col] = cur_col.size();
		final_mat->data[col] = new pairIF[cur_col.size()];

		for(int i = 0; i < cur_col.size(); ++i)
			final_mat->data[col][i] = cur_col[i];
	}
	return final_mat;
}

SMatF* restrict_mat( SMatF* mat, SMatF* pos_mat)
{
	_int cols = pos_mat->nc;	
	SMatF* ans = new SMatF(pos_mat);
	for(_int col=0;col<cols;col++)
	{
		for(_int i=0; i < pos_mat->size[col]; i++)
		{
			_float val = getvalue(pos_mat->data[col][i].first, col, mat);
			ans->data[col][i].second = val;
		}
	}
	ans = eliminate_zeros(ans);
	return ans;
}

SMatF* my_func( SMatF* score_mat, SMatF* Y_token, _float alpha )
{
	cout << alpha <<endl;
	SMatF* ans = new SMatF(Y_token);
	Y_token = make_binary(Y_token); 	
	Y_token = Y_token->transpose();
	cout << Y_token->nr << Y_token->nc <<endl;
	cout << score_mat->nr << score_mat->nc <<endl;
	SMatF* f_score_mat = transform_matf(score_mat,alpha);
	SMatF* g_score_mat = transform_matg(score_mat,alpha);

	SMatF* pos_mat = non_zero_prod_pos_mat(Y_token, score_mat);
	cout << "reached1 "<<endl;	
	SMatF* f_prodmat = Y_token->prod(f_score_mat);
	cout << "reached2 "<<endl;
	SMatF* g_prodmat = Y_token->prod(g_score_mat);
	cout << "reached3 "<<endl;
	f_prodmat = restrict_mat(f_prodmat,pos_mat);
	g_prodmat = restrict_mat(g_prodmat,pos_mat);
	cout << "reached "<<endl;
	for(_int col=0; col<ans->nc; col++)	
	{
		if(f_prodmat->size[col]!=g_prodmat->size[col])
			cout<< f_prodmat->size[col] << "  " << g_prodmat->size[col] <<endl;
	}
	/*
		for(_int i=0;i<ans->size[col];i++)
		{
			if(f_prodmat->data[col][i].first != g_prodmat->data[col][i].first)
			{
				cout << "error here " << f_prodmat->data[col][i].first <<"  "<< g_prodmat->data[col][i].first <<endl;
			}
			//ans->data[col][i].second = f_prodmat->data[col][i].second/g_prodmat->data[col][i].second;
		}
*/		
/*
SMatF* non_zero_prod_pos_mat(SMatF* mat1, SMatF* mat2)
{
	SMatF* one_mat2 = make_binary(mat2);
	SMatF* lbl_prodmat = mat1->prod(one_mat2);
	SMatF* mat1_t = mat1->transpose();
	SMatF* lbl_pos_mat = new SMatF(lbl_prodmat->nr,lbl_prodmat->nc);

	for(_int col=0;col<lbl_pos_mat->nc;col++)
	{
		_int siz=0;
		vector< pair<_int, _float> > vec;//vector<pair<_int, _float>>
		for(_int i=0;i<lbl_prodmat->size[col];i++)
		{
			if(lbl_prodmat->data[col][i].second == mat1_t->size[col])	
			{
				siz++;
				vec.push_back(pairIF(lbl_prodmat->data[col][i].first, 1.0) );				
			}
		}

		lbl_pos_mat->data[col]=new pair<_int,_float>[siz];
		for(_int i=0;i<siz;i++)
			lbl_pos_mat->data[col][i]=vec[i];

		lbl_pos_mat->size[col]=siz;
	}
	return lbl_pos_mat;
}

SMatF* all_non_zero_prod(SMatF* mat1, SMatF* mat2)
{
	SMatF* lbl_pos_mat = non_zero_prod_pos_mat(mat1,mat2);
	SMatF* prodmat = mat1->sparse_prod(mat2,lbl_pos_mat);
	return prodmat;
}
*/

/*
SMatF* my_func( SMatF* score_mat, SMatF* Y_token, _float alpha )
{
	cout << alpha <<endl;
	
	Y_token = make_binary(Y_token); 	
	SMatF* f_score_mat = transform_matf(score_mat,alpha);
	SMatF* g_score_mat = transform_matg(score_mat,alpha);
	cout << "reached here"<<endl;
	_int dps = score_mat -> nc;
	_int labs = Y_token -> nc;

	SMatF* final_mat = new SMatF(labs, dps);

	for(_int dp=0;dp<dps;dp++)
	{
		//pair<_int,_float>* vec_dp = f_score_mat[dp];
		vector<pairIF> cur_col;
		for(_int lab=0;lab<labs;lab++)
		{		
			_int toks = Y_token->size[lab];	
			pair<_int,_float>* vec_lab = Y_token->data[lab];	
			_float prodf = 0, prodg = 0;
			_int f_ind = 0, g_ind=0;
			for(_int tok=0;tok<toks;tok++)
			{
				f_ind = find_ind(vec_lab[tok].first,dp,f_score_mat,f_ind);
				//_float val = getvalue(vec_lab[tok].first,dp,f_score_mat);
				if(f_ind == -1)
				{
					prodf=0;
					break;	
				}
				/*
				if(abs(val) < 1e-10)
				{
					prodf=0;
					break;
				}
				*/
/*
				_float f_val = f_score_mat->data[dp][f_ind].second;
				prodf += vec_lab[tok].second * f_val;

				g_ind = find_ind(vec_lab[tok].first,dp,g_score_mat,g_ind);
				_float g_val = g_score_mat->data[dp][g_ind].second;
				prodg += vec_lab[tok].second * g_val;				
				//prodg += vec_lab[tok].second*getvalue(vec_lab[tok].first,dp,g_score_mat);
			}
			if(abs(prodf) > 1e-10)
				cur_col.push_back(pairIF(lab, prodf / prodg));
		}

		sort(cur_col.begin(), cur_col.end(), comp_pair_by_second_desc<_int,_float>);
		cur_col.resize(min((size_t)100, cur_col.size()));

		final_mat->size[dp] = cur_col.size();
		final_mat->data[dp] = new pairIF[cur_col.size()];

		for(int i = 0; i < cur_col.size(); ++i)
			final_mat->data[dp][i] = cur_col[i];
	}

	SMatF* ans = transform_exp(final_mat);
	return ans;
}
*/


